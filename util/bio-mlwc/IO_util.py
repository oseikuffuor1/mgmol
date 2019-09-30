import sys
import numpy as np
from map_util import nearest_atoms, MLWC

# This module contains methods for reading, writing, and altering files
# in formats that are common within MGmol and for biological problems.

class Atom:
	"""
	Class storing ATOM or HETATM records from a .pdb file.

	param name (string): atom name
	param res_name (string): residue name
	param res_num (int): residue number
	param coords (np.array): length 3 array of coordinates
	param element (string): element name
	"""

	def __init__(self, name, res_name, res_num, coords, element):
		self.name = name
		self.res_name = res_name
		self.res_num = res_num
		self.coords = coords
		self.element = element

	def __repr__(self):
		return "Atom({}, residue={}#{}, element={}, {})".format(self.name, self.res_name, self.res_num, self.element, [round(c, 2) for c in self.coords])


class SSbond:
	"""
	Class storing SSBOND records from a .pdb file.

	param res_num1 (int): residue number of residue involved in SSBOND
	param res_num2 (int): residue number of residue involved in SSBOND
	"""

	def __init__(self, res_num1, res_num2):
		self.res_num1 = res_num1
		self.res_num2 = res_num2


def read_atoms(file_handle):
	"""
	Reads all ATOM and HETATM records from a file into a list of Atom objects

	param file_handle (file): open .pdb file
	return (list(Atom)): list of Atoms
	"""

	lines = file_handle.readlines()

	atoms = []

	for line in lines:
		if 'ATOM' in line[0:6] or 'HETATM' in line[0:6]:
			atom = Atom(
				name=line[12:16].strip(),
				res_name=line[17:21].strip(),
				res_num=int(line[22:26]),
				coords=np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
				element=line[76:78].strip()
				)
			atoms.append(atom)

	return atoms


def read_ssbonds(file_handle):
	"""
	Reads all SSBOND records from a .pdb file into a list of SSbond objects

	param file_handle (file): open .pdb file
	return (list(SSbond)): list of SSbonds
	"""

	lines = file_handle.readlines()

	ssbonds = []

	for line in lines:
		if 'SSBOND' in line[0:6]:
			ssbond = SSbond(
				res_num1=int(line[17:21]),
				res_num2=int(line[31:35])
				)
			ssbonds.append(ssbond)

	return ssbonds


def read_xyz(file_handle):
	"""
	Reads coordinates from a .xyz file into a list numpy arrays

	param file_handle (file): open .xyz file
	return (list(np.array)): list of length 3 coordinate numpy arrays
	"""

	lines = file_handle.readlines()[2:]

	points = []

	for line in lines:
		point = np.array([float(field) for field in line.split()[1:4]])
		points.append(point)

	return points


def read_mlwcs(file_handle):
	"""
        Reads coordinates from a .xyz file into a list of MLWC objects

        param file_handle (file): open .xyz file
        return (list(np.array)): list of MLWC objects
        """

	return [MLWC(coords, None) for coords in read_xyz(file_handle)] 


def mlwc_xyz_string(mlwcs, annotate=False):
	"""
	Generates a string formatted as a .xyz file containing the coordinates of all input MLWCs

	param mlwcs (list(np.array)): list of length 3 coordinate numpy arrays
	param annotate (bool): whether or not to include info on how each MLWC was generated
	return (string): xyz-formatted string containing all coordinates
	"""

	s = "{}\nMLWC\n".format(len(mlwcs))
	for mlwc in mlwcs:
		if annotate:
			s += 'Wa'+' '+str(mlwc.coords[0])+' '+str(mlwc.coords[1])+' '+str(mlwc.coords[2])+"\t\t# "+str(mlwc.info)+'\n'
		else:
			s += 'Wa'+' '+str(mlwc.coords[0])+' '+str(mlwc.coords[1])+' '+str(mlwc.coords[2])+'\n'

	return s


def mlwc_lrs_string(mlwcs):
	"""
	Generates a string formatted as a .lrs file containing the coordinates of all input MLWCs

	param mlwcs (list(np.array)): list of length 3 coordinate numpy arrays
	return (string): in-formatted string containing all coordinates
	"""

	ang_to_bohr = 1.889725989
	s = ""

	for mlwc in mlwcs:
		s += str(ang_to_bohr * mlwc.coords[0])+'\t'+str(ang_to_bohr * mlwc.coords[1])+'\t'+str(ang_to_bohr * mlwc.coords[2])+'\n'

	return s


def clean_pdb(file_handle):
	"""
	Removes multiple conformations and non-essential information from a .pdb file

	param file_handle (file): open .pdb file
	return (string): pdb-formatted string with all non-redundant atoms and disulfide bonds
	"""

	lines = file_handle.readlines()

	pdb = ""

	for line in lines:
		# extract only atoms and ssbonds
		# don't take beta (or other) conformations - it results in duplicate atoms
		if any([i in line[0:6] for i in ['ATOM', 'HETATM', 'SSBOND']]) and line[16] in ["A", " "]:
			pdb += line

	return pdb


def transfer_geom(xyz_file_handle, pdb_file_handle):
	"""
	Takes the ionic positions in a .xyz file and applies them to the atoms in a .pdb file, preserving all other fields.
	This can be used after running mgmol2xyz.py on the output of geometry optimization to get optimized ionic positions.
	Then the MLWCs can be geometrically generated for pdb,
	which now has the same residue information, just with optimized geometry.

	param xyz_file_handle (file): open .xyz file
	param pdb_file_handle (file): open .pdb file
	return (string): pdb-formatted string with all information from .pdb preserved except atomic coordinates taken from .xyz
	"""

	xyz_lines = xyz_file_handle.readlines()[2:]
	pdb_lines = [line for line in pdb_file_handle.readlines() if 'ATOM' in line[0:6] or 'HETATM' in line[0:6]]

	if len(xyz_lines) != len(pdb_lines):
		print "\nERROR: different number of atoms found in .xyz ({}) and .pdb ({}) files\n".format(len(xyz_lines), len(pdb_lines))
		exit()

	pdb = ""

	for i in range(0, len(pdb_lines)):
		xyz_line = xyz_lines[i]
		pdb_line = pdb_lines[i]

		coords = [float(num) for num in xyz_line.split()[1:]]
		pdb_line = pdb_line[:30] + "".join([padded_num(n) for n in coords]) + pdb_line[54:]
		pdb += pdb_line

	return pdb


def padded_num(n):
	"""
	Pads a number to be 8 characters long in order to be inserted in a .pdb file

	param n (num): number to be padded
	return (string): 8-character string containing padded, possibly rounded number
	"""

	s = str(round(n, 3))
	s = " "*(8-len(s)) + s

	return s


if __name__ == "__main__":
	try:
		command = sys.argv[1]

		if command == "clean_pdb":
			in_file_name = sys.argv[2]
			out_file_name = sys.argv[3]
		elif command == "transfer_geom":
			geom_file_name = sys.argv[2]
			in_file_name = sys.argv[3]
			out_file_name = sys.argv[4]
		else:
			raise IndexError

	except IndexError:
		print "\nUSAGE: python {} clean_pdb <input.pdb> <output.pdb>\n".format(sys.argv[0])
		print "USAGE: python {} transfer_geom <geom.xyz> <input.pdb> <output.pdb>\n".format(sys.argv[0])
		exit()

	if command == "clean_pdb":
		with open(in_file_name, "r") as file_handle, open(out_file_name, "w") as out_file:
			out_file.write(clean_pdb(file_handle))
	elif command == "transfer_geom":
		with open(geom_file_name, "r") as xyz_file_handle, open(in_file_name, "r") as pdb_file_handle, open(out_file_name, "w") as out_file:
			out_file.write(transfer_geom(xyz_file_handle, pdb_file_handle))
