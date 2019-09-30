import sys
import numpy as np
from IO_util import read_xyz


def accuracy_stats(mlwcs1, mlwcs2):
	ang_to_bohr = 1.889725989

	max_vec = None
	max_norm = None
	max_indices = None

	min_vec = None
	min_norm = float('inf')
	min_indices = None

	total_dist = 0

	for i in range(0, len(mlwcs1)):
		c1 = mlwcs1[i]
		j = min(range(0, len(mlwcs2)), key=lambda j: np.linalg.norm(c1 - mlwcs2[j]))
		c2 = mlwcs2[j]
		diff = np.linalg.norm(c1 - c2)

		if diff > max_norm:
			max_vec = c1 - c2
			max_norm = diff
			max_indices = (i,j)
		if diff < min_norm:
			min_vec = c1 - c2
			min_norm = diff
			min_indices = (i,j)

		total_dist += diff

	return (
		ang_to_bohr * total_dist/len(mlwcs1),
		ang_to_bohr * min_norm,
		ang_to_bohr * max_norm,
		max_indices,
		ang_to_bohr * max_vec
	)


def print_accuracy_stats(stats, file_name1, file_name2):
	avg_norm, min_norm, max_norm, max_indices, max_vec = stats

	print "\nresults (in Bohr):"
	print "\navg |c_1 - c_2|: {}".format(avg_norm)
	print "min |c_1 - c_2|: {}".format(min_norm)
	print "max |c_1 - c_2|: {}".format(max_norm)
	print "\nlargest difference:"
	print "MLWC {} from {}".format(max_indices[0], file_name1)
	print "MLWC {} from {}".format(max_indices[1], file_name2)
	print "c_1 - c_2 = {}, {}, {}\n".format(max_vec[0], max_vec[1], max_vec[2])


def redundancy_stats(mlwcs):
	ang_to_bohr = 1.889725989

	min_vec = None
	min_norm = float('inf')
	min_indices = None

	for i in range(0, len(mlwcs)):
		for j in range(0, len(mlwcs)):
			if i != j:
				c1 = mlwcs[i]
				c2 = mlwcs[j]
				diff = np.linalg.norm(c1 - c2)

				if diff < min_norm:
					min_vec = c1 - c2
					min_norm = diff
					min_indices = (i,j)

	return ang_to_bohr * min_norm, min_indices, ang_to_bohr * min_vec


def print_redundancy_stats(stats):
	min_norm, min_indices, min_vec = stats

	print "\nresults (in Bohr):"
	print "\nclosest MLWCs:"
	print "MLWC {}".format(min_indices[0])
	print "MLWC {}".format(min_indices[1])
	print "|c_1 - c_2|: {}".format(min_norm)
	print " c_1 - c_2 = {}, {}, {}\n".format(min_vec[0], min_vec[1], min_vec[2])


if __name__ == "__main__":
	try:
		file_name1 = sys.argv[1]
		file_name2 = sys.argv[2]

	except IndexError:
		print "\nUSAGE: python {} <input1.xyz> <input2.xyz>\n".format(sys.argv[0])
		print "if the file names provided are the same, statistics regarding the closest, possibly redundant MLWCs will be printed"
		print "otherwise, statistics regarding the largest difference in the location of a given MLWC between the two files will be printed\n"
		exit()

	with open(file_name1, "r") as f1, open(file_name2, "r") as f2:
		mlwcs1 = read_xyz(f1)
		mlwcs2 = read_xyz(f2)

	if len(mlwcs1) != len(mlwcs2):
		print "\nERROR: different number of MLWCs found in {} ({}) and {} ({})\n".format(file_name1, len(mlwcs1), file_name2, len(mlwcs2))
		exit()

	# if the files are the same, give the closest MLWCs (they might be redundant or linearly dependent)
	if file_name1 == file_name2:
		print_redundancy_stats(redundancy_stats(mlwcs1))

	# otherwise give statistics on the distance of each MLWC in file 1 to the nearest MLWC in file 2
	# the nearest MLWC is assumed to be the same MLWC because order is not guaranteed when the MLWCs are generated differently
	else:
		print_accuracy_stats(accuracy_stats(mlwcs1, mlwcs2))
