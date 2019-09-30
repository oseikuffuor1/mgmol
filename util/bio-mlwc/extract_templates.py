import sys
import os
from map_util import nearest_atoms
from IO_util import read_atoms, read_mlwcs, mlwc_xyz_string

# This script contains the methods necessary to extract template files
# from the output of MGmol's Quench routine.

def extract_residue_pdb(pdb_lines, out_file_name):
	"""
	Extracts .pdb file for the second residue in the input. This is useful when
	Quench has been run on tripeptides or 3-base DNA strands.

	param pdb_lines (list(str)): lines of .pdb file
	param out_file_name (string): name of file to contain single-residue .pdb lines WITHOUT EXTENSION
	"""

	residue_pdb = "HEADER\n"

	for line in pdb_lines:
		if 'ATOM' in line[0:6] or 'HETATM' in line[0:6]:
			atom_name = line[12:16].strip()
			residue_number = float(line[22:26])

			if residue_number == 2:
				residue_pdb += line

	with open("{}.pdb".format(out_file_name), "w") as pdb_file:
		pdb_file.write(residue_pdb)


def extract_residue_xyz(mlwcs, atoms, out_file_name):
	"""
	Extracts .xyz file for MLWCs of the second residue in the input. This is useful when
	Quench has been run on tripeptides or 3-base DNA strands.

	param mlwcs (list(np.array)): list of length 3 numpy arrays of MLWC coordinates
	param atoms (list(Atom)): list of Atom objects
	param out_file_name (string): name of file to contain single-residue MLWCs WITHOUT EXTENSION
	"""

	# extract xyz files for the MLWCs of the middle residue by proximity to atoms
	residue_mlwcs = []

	for mlwc in mlwcs:
		bond_atoms = nearest_atoms(mlwc, atoms)[:2]

		# associate MLWC with nearest residue
		nearest_residue_1 = bond_atoms[0].res_name
		nearest_residue_2 = bond_atoms[1].res_name
		nearest_residue_number = bond_atoms[0].res_num
		if nearest_residue_1 == nearest_residue_2 and nearest_residue_number == 2:
			residue_mlwcs.append(mlwc)

	with open("{}-mlwc.xyz".format(out_file_name), "w") as xyz_file:
		xyz_file.write(mlwc_xyz_string(residue_mlwcs))


def extract_all_residue_templates(in_file_names, out_file_names):
	"""
	Extracts a .pdb file for ionic positions and a .xyz file for MLWCs of the second residue in every input.
	This is useful when Quench has been run on tripeptides or 3-base DNA strands.

	param in_file_names (list(string)): list of filenames containing tripeptide ions and MLWCs WITHOUT EXTENSION
	param out_file_names (list(string)): list of filenames to contain single-residue ions and MLWCs WITHOUT EXTENSION
	"""

	for i in range(0, len(in_file_names)):
		in_file_name = in_file_names[i]
		out_file_name = out_file_names[i]
		tripeptide_atom_file_name = "{}.pdb".format(in_file_name)
		tripeptide_mlwc_file_name = "{}-mlwc.xyz".format(in_file_name)

		with open(tripeptide_atom_file_name, "r") as tripeptide_atom_file, open(tripeptide_mlwc_file_name, "r") as tripeptide_mlwc_file:
			extract_residue_pdb(tripeptide_atom_file.readlines(), out_file_name)
			tripeptide_atom_file.seek(0)
			extract_residue_xyz(
				read_mlwcs(tripeptide_mlwc_file),
				read_atoms(tripeptide_atom_file),
				out_file_name
			)


def extract_peptide_bond_pdb(pdb_lines, out_file_name):
	"""
	Extracts .pdb file for the peptide bond.

	param pdb_lines (list(str)): lines of .pdb file
	param out_file_name (string): name of file to contain peptide bond .pdb WITHOUT EXTENSION
	"""

	peptide_bond_pdb = "HEADER\n"

	for line in pdb_lines:
		if 'ATOM' in line[0:6] or 'HETATM' in line[0:6]:
			atom_name = line[12:16].strip()
			residue_number = float(line[22:26])

			# extract peptide bond template from between first two residues
			if (atom_name == "C" and residue_number == 1) or (atom_name in ["N", "CA"] and residue_number == 2):
				peptide_bond_pdb += line

	with open("{}.pdb".format(out_file_name), "w") as peptide_bond_file:
		peptide_bond_file.write(peptide_bond_pdb)


def extract_peptide_bond_xyz(mlwcs, atoms, out_file_name):
	"""
	Extracts .xyz file for MLWCs of the peptide bond.

	param mlwcs (list(np.array)): list of length 3 numpy arrays of MLWC coordinates
	param atoms (list(Atom)): list of Atom objects
	param out_file_name (string): name of file to contain peptide bond MLWCs WITHOUT EXTENSION
	"""

	bond_mlwcs = []

	for mlwc in mlwcs:
		bond_atoms = nearest_atoms(mlwc, atoms)[:2]
		bond_atom_names = [a.name for a in bond_atoms]

		# extract peptide bond template from between first two residues
		if all([name in ["C", "N"] for name in bond_atom_names]):
			res_num_1 = bond_atoms[0].res_num
			res_num_2 = bond_atoms[1].res_num
			if res_num_1 in [1, 2] and res_num_2 in [1, 2]:
				bond_mlwcs.append(mlwc)

	with open("{}-mlwc.xyz".format(out_file_name), "w") as peptide_bond_file:
		peptide_bond_file.write(mlwc_xyz_string(bond_mlwcs))


def extract_pde_bond_pdb(pdb_lines, out_file_name):
	"""
	Extracts .pdb file for the phosphodiester bond.

	param pdb_lines (list(str)): lines of .pdb file
	param out_file_name (string): name of file to contain phosphodiester bond .pdb WITHOUT EXTENSION
	"""

	# phosphodiester bond pdb
	pde_bond_pdb = "HEADER\n"

	for line in pdb_lines:
		if 'ATOM' in line[0:6] or 'HETATM' in line[0:6]:
			atom_name = line[12:16].strip()
			residue_number = float(line[22:26])

			# extract phosphodiester bond template from between first two residues
			if (atom_name in ["O3*", "C3*"] and residue_number == 1) or (atom_name == "P" and residue_number == 2):
				pde_bond_pdb += line

	with open("{}.pdb".format(out_file_name), "w") as pde_bond_file:
		pde_bond_file.write(pde_bond_pdb)


def extract_pde_bond_xyz(mlwcs, atoms, out_file_name):
	"""
	Extracts .xyz file for MLWCs of the phosphodiester bond.

	param mlwcs (list(np.array)): list of length 3 numpy arrays of MLWC coordinates
	param atoms (list(Atom)): list of Atom objects
	param out_file_name (string): name of file to contain phosphodiester bond MLWCs WITHOUT EXTENSION
	"""

	# phosphodiester bond MLWCs
	bond_mlwcs = []

	for mlwc in mlwcs:
		bond_atoms = nearest_atoms(mlwc, atoms)[:2]
		bond_atom_names = [a.name for a in bond_atoms]

		# extract peptide bond template from between first two residues
		if all([name in ["O3*", "P"] for name in bond_atom_names]):
			res_num_1 = bond_atoms[0].res_num
			res_num_2 = bond_atoms[1].res_num
			if res_num_1 in [1, 2] and res_num_2 in [1, 2]:
				bond_mlwcs.append(mlwc)

	with open("{}-mlwc.xyz".format(out_file_name), "w") as pde_bond_file:
		pde_bond_file.write(mlwc_xyz_string(bond_mlwcs))


if __name__ == "__main__":
	try:
		system_type = sys.argv[1]
	except IndexError:
		print "\nUSAGE: python {} <protein|DNA|RNA>\n".format(sys.argv[0])
		exit()

	if system_type == "protein":
		residues = [
		"ALA", "ARG", "ASN", "ASP",
		"CYS", "GLN", "GLU", "GLY",
		"VAL", "ILE", "LEU", "LYS",
		"MET", "PHE", "PRO", "SER",
		"THR", "TRP", "TYR", "HSP",
		"HSD", "HSE"
		]
		in_file_names = ["protein/template_files/{}/{}-tripeptide".format(r, r) for r in residues]

		with open("protein/template_files/disulfide.pdb", "r") as atom_file, open("protein/template_files/disulfide-mlwc.xyz", "r") as mlwc_file:
			extract_peptide_bond_pdb(atom_file.readlines(), "protein/template_files/peptide_bond/peptide_bond")
			atom_file.seek(0)
			extract_peptide_bond_xyz(
				read_mlwcs(mlwc_file),
				read_atoms(atom_file),
				"protein/template_files/peptide_bond/peptide_bond"
			)

	elif system_type in ["DNA", "RNA"]:
		if system_type == "DNA":
			residues = ["A", "C", "G", "T"]
			in_file_names = ["DNA/template_files/{}/dna_a{}a".format(r, r.lower()) for r in residues]

		elif system_type == "RNA":
			residues = ["A", "C", "G", "U"]
			in_file_names = ["RNA/template_files/{}/rna_a{}a".format(r, r.lower()) for r in residues]

		with open("{}.pdb".format(in_file_names[0]), "r") as atom_file, open("{}-mlwc.xyz".format(in_file_names[0]), "r") as mlwc_file:
			extract_pde_bond_pdb(atom_file.readlines(), "{}/template_files/phosphodiester_bond/phosphodiester_bond".format(system_type))
			atom_file.seek(0)
			extract_pde_bond_xyz(
				read_mlwcs(mlwc_file),
				read_atoms(atom_file),
				"{}/template_files/phosphodiester_bond/phosphodiester_bond".format(system_type)
			)

	else:
		print "\nUSAGE: python {} <protein|DNA|RNA>\n".format(sys.argv[0])
		exit()

	out_file_names = ["{}/template_files/{}/{}".format(system_type, r, r) for r in residues]
	extract_all_residue_templates(in_file_names, out_file_names)
