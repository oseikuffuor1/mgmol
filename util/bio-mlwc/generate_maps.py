import os
import sys
import pickle as pkl
import numpy as np
from copy import deepcopy
from IO_util import read_atoms, read_mlwcs
from map_util import generate_mlwc_map, Map

# This module contains the methods necessary to generate Map objects
# from .xyz and .pdb templates once they have been extracted from MGmol outputs.

def generate_amino_acid_map(mlwcs, atoms):
	"""
	Generates a list of maps that take the set of input atoms to the set of MLWCs
	for a given amino acid.

	Assumes symmetry to create maps for the amino and carboxyl charge groups.

	param mlwcs (list(np.array)): list of length 3 numpy arrays of MLWC coordinates associated with this amino acid
	param atoms (list(Atom)): list of Atom objects associated with this amino acid
	return (list(Map)): list of maps associated with this amino acid
	"""

	# generates collection of maps, one for each MLWC in the residue
	mlwc_maps = [generate_mlwc_map(mlwc, atoms) for mlwc in mlwcs]

	# find maps related to the O atom in the carboxyl group
	carboxyl_maps = [m for m in mlwc_maps if m.atom_names[0] == "O" and m.atom_names[1] == "C"]
	for carboxyl_map in carboxyl_maps:
		# duplicate this map, relying on symmetry, to the other O atom in the carboxyl group
		new_carboxyl_map = deepcopy(carboxyl_map)
		new_carboxyl_map.atom_names[0] = "OXT"
		mlwc_maps.append(new_carboxyl_map)

	# find the map related to the H atom in the amino group
	try:
		amino_map = [m for m in mlwc_maps if m.atom_names[0] == "H"][0]
		for i in range(2, 4):
			# duplicate this map, relying on symmetry, to the H2 and H3 atoms in the amino group
			new_amino_map = deepcopy(amino_map)
			new_amino_map.atom_names[0] = str(i) + "H"
			mlwc_maps.append(new_amino_map)
	except IndexError:
		# PRO has no H-N maps to duplicate
		pass

	return mlwc_maps


def generate_nucleobase_map(mlwcs, atoms):
	"""
	Generates a list of maps that take the set of input atoms to the set of MLWCs
	for a given nucleobase.

	Adds a MLWC directly at the center of the O3* - HO3* and O5* - HO5* bonds, which is otherwise missed.

	param mlwcs (list(np.array)): list of length 3 numpy arrays of MLWC coordinates associated with this nucleobase
	param atoms (list(Atom)): list of Atom objects associated with this nucleobase
	return (list(Map)): list of maps associated with this nucleobase
	"""

	# generates collection of maps, one for each MLWC in the residue
	mlwc_maps = [generate_mlwc_map(mlwc, atoms) for mlwc in mlwcs]

	# add center between O3* and HO3*
	mlwc_maps.append(Map(np.array([1,0,0]), 0.5, ["O3*", "HO3*", "X"]))
	# add center between O5* and HO5*
	mlwc_maps.append(Map(np.array([1,0,0]), 0.5, ["O5*", "HO5*", "X"]))

	return mlwc_maps


def generate_simple_map(mlwcs, atoms):
	"""
	Generates a list of maps that take the set of input atoms to the set of MLWCs.

	There's no funny business here - no assumptions of symmetry or added centers.
	We simply create a list of Maps.

	param mlwcs (list(np.array)): list of length 3 numpy arrays of MLWC coordinates
	param atoms (list(Atom)): list of Atom objects
	return (list(Map)): list of maps
	"""

	# generates collection of maps, one for each MLWC
	return [generate_mlwc_map(mlwc, atoms) for mlwc in mlwcs]


def generate_all_maps(residues, residue_atom_file_names, residue_mlwc_file_names, residue_generating_func,
 				      bonds, bond_atom_file_names, bond_mlwc_file_names, bond_generating_func,
					  map_file_name):
	"""
	Generates a list of maps that take the set of input atoms to the set of MLWCs
	for the set of given residues and bonds, whether they be amino acids or nucleobases.

	This is just a helpful abstraction to generalize the map generation process to proteins, DNA, RNA, etc.

	param residues (list(string)): list of residue names
	param residue_atom_file_names (list(string)): list of file names containing the atoms for each residue
	param residue_mlwc_file_names (list(string)): list of file names containing the MLWCs for each residue
	param residue_generating_func (function): function that takes `mlwcs` and `atoms` arguments and returns a list of maps

	param bonds (list(string)): list of bond names
	param bond_atom_file_names (list(string)): list of file names containing the atoms for each bond
	param bond_mlwc_file_names (list(string)): list of file names containing the MLWCs for each bond
	param bond_generating_func (function): function that takes `mlwcs` and `atoms` arguments and returns a list of maps

	param map_file_name (string): name of .pkl file where maps will be dumped
	return (list(list(Map))): list of residue maps, which themselves are lists of Maps
	"""

	maps = {}

	for i in range(0, len(residues)):
		residue = residues[i]
		atom_file_name = residue_atom_file_names[i]
		mlwc_file_name = residue_mlwc_file_names[i]

		with open(atom_file_name, "r") as atom_file, open(mlwc_file_name, "r") as mlwc_file:
			atoms = read_atoms(atom_file)
			mlwcs = read_mlwcs(mlwc_file)

		maps[residue] = residue_generating_func(mlwcs, atoms)

	for i in range(0, len(bonds)):
		bond = bonds[i]
		atom_file_name = bond_atom_file_names[i]
		mlwc_file_name = bond_mlwc_file_names[i]

		with open(atom_file_name, "r") as atom_file, open(mlwc_file_name, "r") as mlwc_file:
			atoms = read_atoms(atom_file)
			mlwcs = read_mlwcs(mlwc_file)

		maps[bond] = bond_generating_func(mlwcs, atoms)

	with open(map_file_name, "wb") as map_file:
		pkl.dump(maps, map_file)

	return maps


if __name__ == "__main__":
	try:
		system_type = sys.argv[1]
	except IndexError:
		print "\nUSAGE: python {} <protein|DNA|RNA|ligands>\n".format(sys.argv[0])
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
		bonds = ["peptide_bond", "disulfide_bond"]
		residue_generating_func = generate_amino_acid_map

	elif system_type == "DNA":
		residues = ["A", "C", "G", "T"]
		bonds = ["phosphodiester_bond"]
		residue_generating_func = generate_nucleobase_map

	elif system_type == "RNA":
		residues = ["A", "C", "G", "U"]
		bonds = ["phosphodiester_bond"]
		residue_generating_func = generate_nucleobase_map

	elif system_type == "ligands":
		residues = []
		# we put the ligands in the `bonds` variable because simple mappings are generated
		# so we don't need to write a redundant MLWC function.
		
		# automatically generate maps for all ligands in template_files
		# this makes it easier to add new ligands just by adding templates
		bonds = [directory for directory in os.listdir("{}/template_files".format(system_type))]
		residue_generating_func = None

	else:
		print "\nUSAGE: python {} <protein|DNA|RNA|ligands>\n".format(sys.argv[0])
		exit()

	residue_atom_file_names = ["{}/template_files/{}/{}.pdb".format(system_type, r, r) for r in residues]
	residue_mlwc_file_names = ["{}/template_files/{}/{}-mlwc-repaired.xyz".format(system_type, r, r)
							   if os.path.isfile("{}/template_files/{}/{}-mlwc-repaired.xyz".format(system_type, r, r))
	 				   		   else "{}/template_files/{}/{}-mlwc.xyz".format(system_type, r, r) for r in residues]

	bond_atom_file_names = ["{}/template_files/{}/{}.pdb".format(system_type, b, b) for b in bonds]
	bond_mlwc_file_names = ["{}/template_files/{}/{}-mlwc-repaired.xyz".format(system_type, b, b)
							if os.path.isfile("{}/template_files/{}/{}-mlwc-repaired.xyz".format(system_type, b, b))
	 				   		else "{}/template_files/{}/{}-mlwc.xyz".format(system_type, b, b) for b in bonds]

	print generate_all_maps(
		residues, residue_atom_file_names, residue_mlwc_file_names, residue_generating_func,
		bonds, bond_atom_file_names, bond_mlwc_file_names, generate_simple_map,
		"{}/{}_maps.pkl".format(system_type, system_type)
		)
