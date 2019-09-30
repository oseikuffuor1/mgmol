import numpy as np
from heapq import nsmallest

# All the real, abstractable math (the geometric mapping procedure) happens here.
# That's why it has its own module apart from the execution scripts.

class MapException(Exception):
	"""
	Exception to be raised when the mapping procedure encounters a difficulty.
	Provides the information necessary to report the error effectively.
	"""

	def __init__(self, message, atom_name=None, res_num=None, res_name=None):
		super(MapException, self).__init__(message)
		self.atom_name = atom_name
		self.res_num = res_num
		self.res_name = res_name


class Map:
	"""
	A map that takes a set of input atoms to the coordinates of a MLWC.

	The map consists of:
	- a vector from the first atom to the MLWC given in a local basis,
	- a length scaling factor so the mapping scales well with various bond lengths,
	- the three atom names needed to reconstruct the local basis in a new input
	"""

	def __init__(self, mlwc_vector_A, length_factor, atom_names):
		self.mlwc_vector_A = mlwc_vector_A
		self.length_factor = length_factor
		self.atom_names = atom_names

	def __repr__(self):
		return "Map({}, {}, {})".format(self.mlwc_vector_A, self.length_factor, self.atom_names)


class MapInfo:
	"""
	Object used to return map information in order to determine how a MLWC was generated.
	This can be extremely useful for debugging and determining map defects.
	"""

	def __init__(self, atoms_used=None, map_used=None, map_residue=None, map_index=None):
		# atoms
		self.atoms_used = atoms_used
		# actual map object
		self.map_used = map_used
		# residue that map belongs to
		self.map_residue = map_residue
		# index within the residue map of map used
		self.map_index = map_index

	def __repr__(self):
		return "map_atoms={}, atoms_used={}, map={}[{}]".format(
			self.map_used.atom_names if self.map_used is not None else None,
			"["+",".join(["{}:{}#{}".format(a.name, a.res_name, a.res_num) for a in self.atoms_used])+"]" if self.atoms_used is not None else None,
			self.map_residue,
			self.map_index
		)


class MLWC:
	"""
	Stores MLWC with mapping information for easy debugging.
	"""

	def __init__(self, coords, info):
		self.coords = coords
		self.info = info


def orthonormal_basis(b, a):
	"""
	Creates an orthonormal basis matrix with b as the first axis,
	with the other axes oriented with respect to a.

	param b (np.array): length 3 numpy array bond vector
	param a (np.array): length 3 numpy array auxiliary vector
	return (np.array): 3x3 numpy array representing the orthonormal basis as a matrix
	"""

	b = b / np.linalg.norm(b)
	a = a / np.linalg.norm(a)

	u_1 = b
	u_2 = np.cross(a, b)
	u_2 = u_2 / np.linalg.norm(u_2)
	u_3 = np.cross(u_1, u_2)
	u_3 = u_3 / np.linalg.norm(u_3)

	A = np.column_stack((u_1, u_2, u_3))

	return A


def nearest_atoms(mlwc, atoms):
	"""
	Finds the 3 nearest atoms to any MLWC

	param coords (MLWC): MLWC object
	param atoms (list(Atom)): list of atom records
	return (list(Atom)): list of 3 atoms nearest the input MLWC
	"""

	return nsmallest(
		3,
		atoms,
		key=lambda atom: np.linalg.norm(mlwc.coords - atom.coords)
		)


def generate_mlwc_map(mlwc, atoms):
	"""
	Generates a map that takes the three nearest atoms to a MLWC.

	The atoms are named bond_1, bond_2, and auxiliary because that's the way they
	are usually ordered by closeness to the MLWC, although it doesn't
	particularly effect the mapping if that's not the case.

	param mlwc (np.array): MLWC for which the map will be generated
	param atoms (list(Atom)): list of atom records that the map may be derived from
	return (np.array, float, list(string)): map
	"""

	nearest = nearest_atoms(mlwc, atoms)
	atom_names = [atom.name for atom in nearest]

	# at least one atom will exist in the template
	# others might not, in which case we pick arbitrary points to form the basis
	bond_1_coords = nearest[0].coords
	bond_2_coords = nearest[1].coords if len(nearest) > 1 else np.array([1000, 0, 0])
	aux_coords    = nearest[2].coords if len(nearest) > 2 else np.array([0, 1000, 0])

	bond_vector = bond_2_coords - bond_1_coords
	aux_vector  = aux_coords - bond_1_coords
	mlwc_vector = mlwc.coords - bond_1_coords

	A = orthonormal_basis(bond_vector, aux_vector)

	# mlwc_vector in basis A
	mlwc_vector_A = np.transpose(A).dot(mlwc_vector)

	length_factor = np.linalg.norm(mlwc_vector) / np.linalg.norm(bond_vector)

	m = Map(mlwc_vector_A, length_factor, atom_names)

	return m


def calculate_mlwc(mlwc_map, atoms, equiv_names):
	"""
	Calculates the coordinates of the MLWC using the given map.
	This is done by attempting to find the atoms with names that are equivalent
	to those found in the map. If they cannot be found, a MapException is raised
	to let the user know that they may be missing some MLWCs, possibly due to atom naming conventions.

	param mlwc_map (np.array, float, list(string)): map
	param atoms (list(Atom)): list of atom records that the map look through to find the necessary atoms
	return (np.array): length 3 numpy array of MLWC coordinates
	"""

	try:
		bond_1_atom = next(atom for atom in atoms if equiv_names(atom.name, mlwc_map.atom_names[0]))
		bond_1_coords = bond_1_atom.coords
	except StopIteration:
		raise MapException(
			"WARNING - atom {} not found in residue {} number {}".format(mlwc_map.atom_names[0], atoms[0].res_name, atoms[0].res_num),
			atom_name=mlwc_map.atom_names[0], res_num=atoms[0].res_num, res_name=atoms[0].res_name
		)

	# if the map contains a secondary atom name, use it.
	# this allows the MLWC's location to be preserved within cylindrical symmetry around the bond
	if len(mlwc_map.atom_names) > 1:
		try:
			bond_2_atom = next(atom for atom in reversed(atoms) if equiv_names(atom.name, mlwc_map.atom_names[1]))
			bond_2_coords = bond_2_atom.coords
		except StopIteration:
			raise MapException(
				"WARNING - atom {} not found in residue {} number {}".format(mlwc_map.atom_names[1], atoms[0].res_name, atoms[0].res_num),
				atom_name=mlwc_map.atom_names[1], res_num=atoms[0].res_num, res_name=atoms[0].res_name
			)
	# if this atom doesn't exist, then it didn't exist in the template so we throw away this symmetry
	else:
		bond_2_atom = None
		bond_2_coords = np.array([1000, 0, 0])

	# if the map contains a tertiary atom name, use it.
	# this allows the MLWC's location to be preserved within the planar structure of the nearby bonds
	if len(mlwc_map.atom_names) > 2:
		try:
			# this `reversed` makes sure different SG atoms are chosen as bond_1 and aux in disulfide bonds
			aux_atom = next(atom for atom in reversed(atoms) if equiv_names(atom.name, mlwc_map.atom_names[2]))
			aux_coords = aux_atom.coords
		except StopIteration:
			# If auxiliary atom is not found, use the origin.
			# This arbitrary choice will discard some rotational information,
			# but will not change the general arrangement of the centers around the bond.
			# It's not worth raising an error and failing to map the MLWC to preserve the aforementioned rotational information.
			aux_atom = None
			aux_coords = np.array([0, 0, 0])
	# if this atom doesn't exist, then it didn't exist in the template so we throw away this symmetry
	else:
		aux_atom = None
		aux_coords = np.array([0, 1000, 0])

	bond_vector = bond_2_coords - bond_1_coords
	aux_vector  = aux_coords - bond_1_coords

	A_i = orthonormal_basis(bond_vector, aux_vector)
	mlwc_vector = A_i.dot(mlwc_map.mlwc_vector_A)
	mlwc_vector = mlwc_map.length_factor * np.linalg.norm(bond_vector) / np.linalg.norm(mlwc_vector) * mlwc_vector

	mlwc_coords = bond_1_coords + mlwc_vector

	return MLWC(
		mlwc_coords,
		MapInfo(
			atoms_used=[atom for atom in [bond_1_atom, bond_2_atom, aux_atom] if atom is not None],
			map_used=mlwc_map
		)
	)
