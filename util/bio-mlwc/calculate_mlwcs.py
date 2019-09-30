import sys
import os
import argparse
import numpy as np
import pickle as pkl
from heapq import nlargest
from copy import deepcopy
from random import sample
from IO_util import read_atoms, read_ssbonds, mlwc_xyz_string, mlwc_lrs_string
from map_util import MapException, MapInfo, MLWC, calculate_mlwc


def calculate_protein_mlwcs(residue_maps, atoms, ssbonds):
	"""
	Calculates all MLWCs for an entire protein given the atoms and ssbonds in the protein.

	param residue_maps (dict): dictionary where keys are residue names and values are residue maps
	param atoms (list(Atom)): list of Atoms in the protein
	param ssbonds (list(SSbond)): list of SSbonds in the protein
	return (list(np.array)): list of length 3 np.arrays giving the coordinates for each MLWC in the protein
	"""

	mlwcs = []

	residue_atoms = []
	prev_residue_number = atoms[0].res_num
	prev_residue = atoms[0].res_name

	prev_C = None
	prev_N = None

	for i in range(0, len(atoms)):
		atom = atoms[i]
		residue = atom.res_name if atom.res_name != "HIS" else "HSP"

		# calculate MLWCs for previous residue
		if atom.res_num != prev_residue_number:
			mlwcs += residue_mlwcs(residue_maps, prev_residue, prev_residue_number, residue_atoms, calculate_amino_acid_mlwcs)
			residue_atoms = []

		residue_atoms.append(atom)
		prev_residue = residue
		prev_residue_number = atom.res_num

		if i == len(atoms)-1:
			mlwcs += residue_mlwcs(residue_maps, prev_residue, prev_residue_number, residue_atoms, calculate_amino_acid_mlwcs)

		# calculate MLWCs for peptide bond
		if equiv_names(atom.name, "C"):
			prev_C = atom
		elif equiv_names(atom.name, "N"):
			prev_N = atom
		# if this is not the first residue in the system AND its not the first residue in a duplicated chain, which would lead
		# to a res_num that is lower than the previous unconnected chain. This avoids mapping centers for non-existent peptide
		# bonds between chains
		elif equiv_names(atom.name, "CA") and prev_C is not None and prev_N is not None and atom.res_num > prev_C.res_num:
			mlwcs += residue_mlwcs(residue_maps, "peptide_bond", prev_residue_number, [atom, prev_C, prev_N], calculate_amino_acid_mlwcs)

	# calculate MLWCs for disulfide bonds
	for atom_pair in ssbond_atoms(atoms, ssbonds):
		# they are sorted so that the SG atoms are read in the correct order for the map
		atom_pair = sorted(atom_pair, key=lambda a: a.res_num)
		ssbond_res_nums = [a.res_num for a in atom_pair]
		# CB atoms from the residues in the SSBOND
		# they are needed as the auxiliary atoms in the map
		CB_atoms = [atom for atom in atoms if equiv_names(atom.name, "CB") and atom.res_num in ssbond_res_nums]
		# they are sorted so that the SG atoms are read in the correct order for the map
		CB_atoms = sorted(CB_atoms, key=lambda a: a.res_num)
		# add MLWCs near SG atoms
		mlwcs += residue_mlwcs(residue_maps, "disulfide_bond", ssbond_res_nums[0], atom_pair + [CB_atoms[0]], calculate_amino_acid_mlwcs)
		mlwcs += residue_mlwcs(residue_maps, "disulfide_bond", ssbond_res_nums[1], list(reversed(atom_pair)) + [CB_atoms[1]], calculate_amino_acid_mlwcs)
		# add MLWC at center of bond
		mlwcs.append(
			MLWC(
				coords=(atom_pair[0].coords + atom_pair[1].coords) / 2,
				info=MapInfo(
					atoms_used=[a for a in atom_pair],
					map_residue="disulfide_bond({},{})".format(ssbond_res_nums[0], ssbond_res_nums[1])
				)
			)
		)

	return mlwcs


def calculate_amino_acid_mlwcs(residue_map, atoms):
	"""
	Calculates all MLWCs for an amino acid given the atoms in the amino acid.

	param residue_map (list(Map)): residue map composed of a list of MLWC maps
	param atoms (list(Atom)): list of Atoms found in the amino acid at hand
	return (list(np.array)): list of length 3 np.arrays giving the coordinates for each MLWC in the amino acid
	"""

	mlwcs = []
	errors = []

	atom_names = [a.name for a in atoms]
	# remedy situations where HA2 and HA3 exist instead of HA1 and HA2
	# this is actually common and must be dealt with
	for i in range(0, len(atoms)):
		atom = atoms[i]
		if all([
			atom.name[0] == "H",
			len(atom.name) > 1 and atom.name[1] in ["A", "B", "G", "D", "E"],
			atom.name[:-1] + "1" not in atom_names,
			atom.name[-1].isdigit() and atom.name[:-1] + str(int(atom.name[-1])+1) not in atom_names
		]):
			atoms[i].name = atom.name[:-1] + "1"

	for j in range(0, len(residue_map)):
		mlwc_map = residue_map[j]
		try:
			mlwc = calculate_mlwc(mlwc_map, atoms, equiv_names)
			mlwc.info.map_index = j
			mlwcs.append(mlwc)
		except MapException as e:
			errors.append(e)

	# consolidate errors to only print those with unique messages
	unique_msgs = set()
	errors = [e for e in errors if e.message not in unique_msgs and not unique_msgs.add(e.message)]

	for e in errors:
		# suppress error messages from charge group maps, which wont be applied
		# to residues that are not the first or last in a sequence
		if not e.atom_name in ["OXT", "2H", "3H"]:
			print e.message

	return mlwcs


def ssbond_atoms(atoms, ssbonds):
	"""
	Returns pairs of Atom objects associated with each SSbond in a protein.

	param atoms (list(Atom)): list of Atoms in the protein
	param ssbonds (list(SSbond)): list of SSbonds in the protein
	return (list(list(Atom))): list of Atom pairs, one pair for each SSbond
	"""

	bond_atoms = []

	for ssbond in ssbonds:
		try:
			atom_1 = next(atom for atom in atoms if equiv_names(atom.name, "SG") and atom.res_num == ssbond.res_num1)
		except StopIteration:
			print "WARNING - disulfide bond atom not found in residue {}".format(ssbond.res_num1)
			# ignore SSBOND - they can't say we didn't warn 'em
			continue

		try:
			atom_2 = next(atom for atom in atoms if equiv_names(atom.name, "SG") and atom.res_num == ssbond.res_num2)
		except StopIteration:
			print "WARNING - disulfide bond atom not found in residue {}".format(ssbond.res_num2)
			# ignore SSBOND - they can't say we didn't warn 'em
			continue

		bond_atoms.append([atom_1, atom_2])

	return bond_atoms


def calculate_nucleic_acid_mlwcs(residue_maps, atoms, ssbonds=None):
	"""
	Calculates all MLWCs for an entire DNA or RNA model given its atoms.

	param residue_maps (dict): dictionary where keys are residue names and values are residue maps
	param atoms (list(Atom)): list of Atoms found in the model
	param ssbonds (None): empty parameter for compatibility with protein MLWC calculation function
	return (list(np.array)): list of length 3 np.arrays giving the coordinates for each MLWC in the nucleic acid
	"""

	mlwcs = []

	all_residues = ["A", "T", "G", "C", "U"]
	num_residues = atoms[-1].res_num

	residue_atoms = []
	prev_residue_number = atoms[0].res_num
	prev_residue = atoms[0].res_name

	prev_C3 = None
	prev_O3 = None

	for i in range(0, len(atoms)):
		atom = atoms[i]
		residue = atom.res_name
		try:
			# DG == G
			residue = all_residues[[r in residue for r in all_residues].index(True)]
		except ValueError:
			print "WARNING - unknown residue {} number {}".format(residue, atom.res_num)
			continue

		# calculate MLWCs for previous residue
		if atom.res_num != prev_residue_number:
			mlwcs += residue_mlwcs(residue_maps, prev_residue, prev_residue_number, residue_atoms, calculate_nucleobase_mlwcs, num_residues=num_residues)
			residue_atoms = []

		residue_atoms.append(atom)
		prev_residue = residue
		prev_residue_number = atom.res_num

		if i == len(atoms)-1:
			mlwcs += residue_mlwcs(residue_maps, prev_residue, prev_residue_number, residue_atoms, calculate_nucleobase_mlwcs, num_residues=num_residues)

		# calculate MLWCs for phosphodiester bond
		if equiv_names(atom.name, "C3*"):
			prev_C3 = atom
		elif equiv_names(atom.name, "O3*"):
			prev_O3 = atom
		# if this is not the first residue in the system AND its not the first residue in a duplicated chain, which would lead
		# to a res_num that is lower than the previous unconnected chain. This avoids mapping centers for non-existent phosphodiester
		# bonds between chains
		elif equiv_names(atom.name, "P") and prev_C3 is not None and prev_O3 is not None and atom.res_num > prev_C3.res_num:
			mlwcs += residue_mlwcs(residue_maps, "phosphodiester_bond", prev_residue_number, [atom, prev_C3, prev_O3], calculate_nucleobase_mlwcs, num_residues=num_residues)

	return mlwcs


def calculate_nucleobase_mlwcs(residue_map, atoms, num_residues):
	"""
	Calculates all MLWCs for a nucleobase given the atoms in the nucleobase.

	param residue_map (list(Map)): residue map composed of a list of MLWC maps
	param atoms (list(Atom)): list of Atoms found in the nucleobase at hand
	param num_residues (int): number of residues in the model - we need this in order to determine when to suppress errors
	return (list(np.array)): list of length 3 np.arrays giving the coordinates for each MLWC in the nucleobase
	"""

	mlwcs = []
	errors = []

	for j in range(0, len(residue_map)):
		mlwc_map = residue_map[j]
		try:
			mlwc = calculate_mlwc(mlwc_map, atoms, equiv_names)
			mlwc.info.map_index = j
			mlwcs.append(mlwc)
		except MapException as e:
			errors.append(e)

	# consolidate errors to only print those with unique messages
	unique_msgs = set()
	errors = [e for e in errors if e.message not in unique_msgs and not unique_msgs.add(e.message)]

	for e in errors:
		# suppress error messages from end group maps, which wont be applied
		# to residues that are not the first in either strand (residue 1 and n/2+1)
		if e.atom_name not in ["HO3*", "HO5*"] and not (e.res_num in [atoms[0].res_num, (num_residues+atoms[0].res_num)/2+1] and any([equiv_names(e.atom_name, a) for a in ["P", "O1P", "O2P"]])):
			print e.message

	return mlwcs


def calculate_ligands_mlwcs(residue_maps, atoms):
	"""
	Calculates all MLWCs for all ligands in a system given their atoms.

	param residue_maps (dict): dictionary where keys are residue names and values are residue maps
	param atoms (list(Atom)): list of Atoms in all ligands
	return (list(np.array)): list of length 3 np.arrays giving the coordinates for each MLWC in the ligands
	"""

	mlwcs = []

	residue_atoms = []
	prev_residue_number = atoms[0].res_num
	prev_residue = atoms[0].res_name

	global ligand_res_dict

	for i in range(0, len(atoms)):
		atom = atoms[i]
		residue = ligand_res_name(atom.res_name) if ligand_res_name(atom.res_name) is not None else atom.res_name

		# calculate MLWCs for previous residue
		if atom.res_num != prev_residue_number:
			mlwcs += residue_mlwcs(residue_maps, prev_residue, prev_residue_number, residue_atoms, calculate_ligand_mlwcs)
			residue_atoms = []

		residue_atoms.append(atom)
		prev_residue = residue
		prev_residue_number = atom.res_num

		if i == len(atoms)-1:
			mlwcs += residue_mlwcs(residue_maps, prev_residue, prev_residue_number, residue_atoms, calculate_ligand_mlwcs)

	return mlwcs


def calculate_ligand_mlwcs(residue_map, atoms):
	"""
	Calculates all MLWCs for a ligand given the atoms in the ligand. Done simply with no bonds or irregularities.

	param residue_map (list(Map)): residue map composed of a list of MLWC maps
	param atoms (list(Atom)): list of Atoms found in the ligand at hand
	param ssbonds (None): empty parameter for compatibility with protein MLWC calculation function
	return (list(np.array)): list of length 3 np.arrays giving the coordinates for each MLWC in the ligand
	"""

	mlwcs = []
	errors = []

	for j in range(0, len(residue_map)):
		mlwc_map = residue_map[j]
		try:
			mlwc = calculate_mlwc(mlwc_map, atoms, equiv_names)
			mlwc.info.map_index = j
			mlwcs.append(mlwc)
		except MapException as e:
			errors.append(e)

	# consolidate errors to only print those with unique messages
	unique_msgs = set()
	errors = [e for e in errors if e.message not in unique_msgs and not unique_msgs.add(e.message)]

	for e in errors:
		print e.message

	return mlwcs


def residue_mlwcs(residue_maps, residue, residue_number, residue_atoms, residue_calc_func, num_residues=None):
	"""
	Generates the MLWCs associated with a specific residue given its atoms, maps, and calculating function.
	Assigns the map_residue for the MLWC.

	param residue_maps (dict): dictionary of string residue names to lists of residue maps
	param residue (str): residue name, key for residue_maps dict
	param residue_number (int): residue number for error reporting purposes
	param residue_atoms (list(Atom)): list of atoms of specific residue
	param residue_calc_func (function): method used to calculate residue MLWCs
	param num_residues (int): number of residues in system, neede for nucleic acids
	return (int): number of valence electrons
	"""

	try:
		if num_residues is not None:
			mlwcs = residue_calc_func(residue_maps[residue], residue_atoms, num_residues=num_residues)
		else:
			mlwcs = residue_calc_func(residue_maps[residue], residue_atoms)

		for mlwc in mlwcs:
			mlwc.info.map_residue = residue

		return mlwcs

	except KeyError:
		print "WARNING - unknown residue {} number {}, no MLWCs mapped".format(residue, residue_number)
		return []


def num_valence(atoms):
	"""
	Determines the number of valence electrons in the input. This allows us to determine how many
	doubly occupied orbitals - how many MLWCs - we need.

	param atoms (list(Atoms)): atoms in the input
	return (int): number of valence electrons
	"""

	n = 0

	for atom in atoms:
		if atom.element in ["H"]:
			n += 1
		elif atom.element in ["C"]:
			n += 4
		elif atom.element in ["N", "P"]:
			n += 5
		elif atom.element in ["O", "S"]:
			n += 6
		elif atom.element in ["K"]:
			n += 9
		# some .pdb files don't have element designations, so try to guess from the atom.name
		elif any([l == atom.name[0] for l in ["H"]]):
			n += 1
		elif any([l == atom.name[0] for l in ["C"]]):
			n += 4
		elif any([l == atom.name[0] for l in ["N", "P"]]):
			n += 5
		elif any([l == atom.name[0] for l in ["O", "S"]]):
			n += 6
		elif any([l == atom.name[0] for l in ["K"]]):
			n += 9
		else:
			print "WARNING - no valence known for element {} atom name {} residue {} number {}".format(atom.element, atom.name, atom.res_name, atom.res_num)
			print "        - this atom is being skipped in calculating the number of orbitals, which may result in too few orbitals"

	return n


def residue_sequence(atoms):
	"""
	Extracts the sequence of residues in the input. We need this to estimate the charge, as we know
	which residues are charged.

	param atoms (list(Atoms)): atoms in the input
	return (list(string)): list of residue names
	"""

	prev_residue = atoms[0].res_name
	prev_residue_number = atoms[0].res_num
	seq = [prev_residue]

	for atom in atoms:
		if atom.res_num != prev_residue_number and ligand_res_name(atom.res_name) is None:
			seq.append(atom.res_name)
		prev_residue_number = atom.res_num

	return seq


def system_charge(atoms, system_type):
	"""
	Determines the charge of a given system.

	param atoms (list(Atoms)): atoms in the input
	param system_type (protein|DNA|RNA): system type
	return (int): estimated charge of the system
	"""

	seq = residue_sequence(atoms)

	if system_type == "protein":
		# add the number of bases, subtract the acids
		charge = sum(res_name in ["LYS", "ARG"] for res_name in seq)
		charge -= sum(res_name in ["ASP", "GLU"] for res_name in seq)

	elif system_type in ["DNA", "RNA"]:
		# 3' ends are uncharged, every other nucleobase is negatively charged
		charge = 2 - len(seq)

	return charge


def fill_mlwcs(atoms, geom_mlwcs, system_type, addition_type, removal_type, num_mlwcs):
	"""
	Generates or removes any MLWCs necessary to obtain the number calculated by valence.

	param atoms (list(Atoms)): atoms in the input
	param geom_mlwcs (list(np.array)): list of MLWCs generated by the geometric mapping approach
	param system_type (protein|DNA|RNA): system type
	param addition_type (far|None): method of adding necessary MLWCs
	param removal_type (near|far|None): method of removing access MLWCs
	param num_mlwcs (int): total number of MLWCs necessary
	return (list(np.array), int, int): list of geometric MLWCs augmented with any additional neceassry MLWCs
	"""

	mlwcs = deepcopy(geom_mlwcs)

	# the correct number of MLWCs were generated
	if num_mlwcs == len(mlwcs):
		return mlwcs

	# too few were generated, add MLWCs if flagged
	elif num_mlwcs > len(mlwcs):
		if addition_type == "far":
			# find the atom farthest away from a MLWC and add a MLWC there
			farthest_atoms = nlargest(
				num_mlwcs-len(mlwcs),
				atoms,
				key=lambda a: min([np.linalg.norm(a.coords - mlwc.coords) for mlwc in mlwcs]) if len(mlwcs) > 0 else 0
				)
			mlwcs += [MLWC(a.coords, MapInfo(atoms_used=a)) for a in farthest_atoms]

		elif addition_type == "random":
			# add MLWCs at atom centers sampled at random from all atoms
			mlwcs += [MLWC(a.coords, MapInfo(atoms_used=a)) for a in sample(atoms, num_mlwcs-len(mlwcs))]

	# too many were generated, subtract MLWCs
	else:
		if removal_type == "near":
			# remove MLWCs close to other MLWCs as they are likely to be redundant
			for i in range(0, len(mlwcs)-num_mlwcs):
				index_to_remove = min(
					range(0, len(mlwcs)),
					key=lambda i: min([np.linalg.norm(mlwcs[i].coords - mlwc.coords) for mlwc in mlwcs])
					)
				del mlwcs[index_to_remove]

		elif removal_type == "far":
			# remove MLWCs farthest from atoms as they are likely to be non-physical
			indices_to_remove = nlargest(
				len(mlwcs)-num_mlwcs,
				range(0, len(mlwcs)),
				key=lambda i: min([np.linalg.norm(a.coords - mlwcs[i].coords) for a in atoms])
				)
			for i in sorted(indices_to_remove, reverse=True):
				del mlwcs[i]

		elif removal_type == "random":
			# remove MLWCs at random sampled from all MLWCs
			indices_to_remove = sample(range(0, len(mlwcs)), len(mlwcs)-num_mlwcs)
			for i in sorted(indices_to_remove, reverse=True):
				del mlwcs[i]

	return mlwcs


def equiv_names(n1, n2):
	"""
	Determines if two atom names are equivalent. The naming conventions for atoms in .pdb files
	is highly non-uniform, so this is needed in order to apply the maps to atoms that are equivalent,
	despite having non-identical names.

	param n1 (string): name of atom 1
	param n2 (string): name of atom 2
	return (bool): whether the atom names are equivalent
	"""

	return _equiv_names(n1, n2) or _equiv_names(n2, n1)
def _equiv_names(n1, n2):
	### general equivalences
	if n1 == n2:
		return True
	if n1[-1].isdigit():
		if n1[-1] == "1":
			if n1[:-1] == n2:
				# HG1 == HG
				return True
			if n1[1:-1] + n1[0] == n2:
				# 3HD1 == HD3
				return True
		if n1[-1] + n1[:-1] == n2:
			# HG23 == 3HG2
			return True
	if n1[0].isdigit():
		if n1[0] == "1" and n1[1:] == n2:
			# 1HG2 == HG2
			return True

	### protein equivalences
	if all(n in ["HN", "H", "1H", "HT1"] for n in [n1, n2]):
		return True
	if all(n in ["2H", "HT2"] for n in [n1, n2]):
		return True
	if all(n in ["3H", "HT3"] for n in [n1, n2]):
		return True
	if all(n in ["O", "OT1"] for n in [n1, n2]):
		return True
	if all(n in ["OXT", "OT2"] for n in [n1, n2]):
		return True

	### nucleic equivalences
	if n1[-1] == "*" and n2[-1] == "'" and equiv_names(n1[:-1], n2[:-1]):
		# 1H5* == H5'
		return True
	if len(n1) >= 2 and n1[-2:] == "''" and n2 == "2{}*".format(n1[:-2]):
		# 2H2* == H2''
		return True
	if all(n in ["O1P", "OP1"] for n in [n1, n2]):
		return True
	if all(n in ["O2P", "OP2"] for n in [n1, n2]):
		return True

	### ligand equivalences
	# water
	if all(n in ["1H", "1HH", "H1"] for n in [n1, n2]):
		return True
	if all(n in ["2H", "2HH", "H2"] for n in [n1, n2]):
		return True
	if all(n in ["OH", "OH2"] for n in [n1, n2]):
		return True
	# potassium ions
	if all(n in ["K", "POT"] for n in [n1, n2]):
		return True

	return False


def ligand_res_name(res_name):
	"""
	Gives the cannonical ligand name (the one stored with the residue map) for any ligand residue name.
	If no equivalent residue name exists, returns None. Can thus be used to determine if an atom is
	a ligand by its res_name.

	param res_name (string): name of residue
	return (string|None): cannonical ligand name or None
	"""

	# dictionary of the cannonical name of each ligand with all possible equivalent residue names
	ligand_res_dict = {
		"HOH": ["HOH", "TIP", "TIP3", "OSP3", "SWM4"],
		"K": ["K", "POT"]
		}

	# check to see if res_name is equivalent to ligand in dictionary
	for ligand, names in ligand_res_dict.iteritems():
		if res_name in names:
			return ligand

	# check to see if res_name matches exactly the name of a ligand in the template directory
	for ligand in os.listdir("ligands/template_files"):
		if res_name == ligand:
			return ligand

	return None


if __name__ == "__main__":
	parser = argparse.ArgumentParser( description='''A utility used to geometrically generate approximate MLWCs for modular biological systems.
			If -a and -r are not specified, no MLWCs will be added or removed.
			MGmol will then add or remove MLWCs sequentially on input, resulting in sections of the system missing MLWCs or with many
			redundant MLWCs.''')
	parser.add_argument('system_type', choices=['protein', 'DNA', 'RNA'], help="system type")
	parser.add_argument('input_file', help="must be a .pdb file")
	parser.add_argument('output_file', help="must be a .xyz or .lrs file")
	parser.add_argument('-c', '--charge', type=int, help="system charge (inferred from residue sequence by default)")
	parser.add_argument('-n', '--num-mlwcs', type=int, help="number of MLWCs to generate (inferred from valence by default)")
	parser.add_argument('-a', '--add', choices=['far', 'random'], help="add MLWCs at the centers of the atoms farthest from other MLWCs, OR at random")
	parser.add_argument('-r', '--remove', choices=['near', 'far', 'random'], help="remove MLWCs nearest other MLWCs (they could be redundant), OR farthest from atoms (they could be unphysical), OR at random")
	parser.add_argument('-A', '--annotate', action='store_true', help="include comments in output file detailing the map and atoms used to generate each MLWC (works only for .xyz files)")

	# parse command line arguments
	args = parser.parse_args()

	system_type = args.system_type
	in_file_name = args.input_file
	out_file_name = args.output_file
	if out_file_name.split(".")[-1] not in ["xyz", "lrs"]:
		parser.print_usage()
		print "calculate_mlwcs.py: error: argument output_file: please use an output with file type `.xyz` or `.lrs`"
		exit()

	if system_type == "protein":
		with open(in_file_name, "r") as infile:
			ssbonds = read_ssbonds(infile)
			infile.seek(0)
		mlwc_calc_func = calculate_protein_mlwcs

	elif system_type in ["DNA", "RNA"]:
		ssbonds = None
		mlwc_calc_func = calculate_nucleic_acid_mlwcs

	else:
		"INVALID SYSTEM TYPE, see -h"
		exit()

	# deserialize system maps for use
	print "\nDeserializing map files {} and ligands/ligands_maps.pkl".format("{}/{}_maps.pkl".format(system_type, system_type))
	try:
		with open("{}/{}_maps.pkl".format(system_type, system_type), "rb") as m_file:
			maps = pkl.load(m_file)
	except IOError:
		print "\nERROR: map file {} does not exist - see README\n".format("{}/{}_maps.pkl".format(system_type, system_type))
		exit()

	# deserialize ligand maps for use
	try:
		with open("{}/{}_maps.pkl".format("ligands", "ligands"), "rb") as m_file:
			ligand_maps = pkl.load(m_file)
	except IOError:
		print "\nERROR: map file {} does not exist - see README\n".format("{}/{}_maps.pkl".format(system_type, system_type))
		exit()

	print "\nReading atoms from file {}".format(in_file_name)
	with open(in_file_name, "r") as infile:
		atoms = read_atoms(infile)
	print "  {} atoms read".format(len(atoms))

	print "\nSeparating out ligand atoms by residue name"
	ligand_atoms = []
	system_atoms = []
	for atom in atoms:
		if ligand_res_name(atom.res_name) is None:
			system_atoms.append(atom)
		else:
			ligand_atoms.append(atom)
	print "  {} ligand atoms found, leaving {} atoms in the biological system".format(len(ligand_atoms), len(system_atoms))

	mlwcs = []

	print "\nCalculating MLWCs for biological system"
	if len(system_atoms) > 0:
		mlwcs += mlwc_calc_func(maps, system_atoms, ssbonds)
	print "  {} MLWCs generated by geometric mapping for the biological system".format(len(mlwcs))

	print "\nCalculating MLWCs for ligands"
	ligand_mlwcs = []
	if len(ligand_atoms) > 0:
		ligand_mlwcs = calculate_ligands_mlwcs(ligand_maps, ligand_atoms)
		mlwcs += ligand_mlwcs
	print "  {} MLWCs generated by geometric mapping for ligands".format(len(ligand_mlwcs))

	print "\n{} MLWCs generated in total by geometric mapping".format(len(mlwcs))

	print "\nAdding atom-centered MLWCs or removing geometrically-generated MLWCs if necessary"
	# take system charge from user if given, otherwise infer from the resiue sequence
	if args.charge is not None:
		charge_method = "given by user"
		system_charge = args.charge
	else:
		charge_method = "inferred from residue sequence"
		system_charge = system_charge(system_atoms + ligand_atoms, system_type)
	# take number of MLWCs from user if given, otherwise use (n+1)/2 where n is the number of valence electrons
	if args.num_mlwcs is not None:
		num_mlwcs_method = "given by user"
		num_mlwcs = args.num_mlwcs
	else:
		num_mlwcs_method = "inferred from system valence"
		num_mlwcs = (num_valence(system_atoms + ligand_atoms) - system_charge + 1) / 2
	print "  system charge: {} ({})".format(system_charge, charge_method)
	print "  MLWCs needed:  {} ({})".format(num_mlwcs, num_mlwcs_method)
	print "  MLWC addition method: {}".format(args.add)
	print "  MLWC removal method:  {}".format(args.remove)
	filled_mlwcs = fill_mlwcs(system_atoms + ligand_atoms, mlwcs, system_type, args.add, args.remove, num_mlwcs)


	if len(mlwcs) > num_mlwcs:
		print "\n{} MLWCs were automatically removed".format(len(mlwcs) - len(filled_mlwcs))
	elif len(mlwcs) < num_mlwcs:
		print "\n{} MLWCs were automatically added".format(len(filled_mlwcs) - len(mlwcs))
	else:
		print "\n0 MLWCs were automatically added or removed"

	if len(filled_mlwcs) > num_mlwcs:
		print "\nThere are {} too many MLWCs in the output. Specify a removal method with -r to generate the correct number of MLWCs.".format(
			len(filled_mlwcs) - num_mlwcs
			)
	elif len(filled_mlwcs) < num_mlwcs:
		print "\nThere are {} too few MLWCs in the output. Specify an addition method with -a to generate the correct number of MLWCs.".format(
			num_mlwcs - len(filled_mlwcs)
			)
	else:
		print "\nThere are the correct number of MLWCs ({}) in the output".format(num_mlwcs)

	mlwcs = filled_mlwcs

	print "\nWriting output to {}\n".format(out_file_name)
	with open(out_file_name, "w") as out_file:
		if out_file_name.split(".")[-1] == "xyz":
			out_file.write(mlwc_xyz_string(mlwcs, annotate=args.annotate))
		elif out_file_name.split(".")[-1] == "lrs":
			out_file.write(mlwc_lrs_string(mlwcs))
