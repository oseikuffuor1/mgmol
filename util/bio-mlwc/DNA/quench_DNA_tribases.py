import os
import sys

# This script runs quench on each tribase using the same .cfg file
# run it from the bio-mlwc/DNA directory
# make sure the bio-mlwc/scratch directory exists so the script can move there to run the computations

if __name__ == "__main__":
	bases = [
	"A", "C", "G", "T"
	]

	file_names = ["../DNA/template_files/{}/dna_a{}a".format(b, b.lower()) for b in bases]

	maindir = "../../.."

	try:
		config_file_name = sys.argv[1]
		if not os.path.isfile(config_file_name):
			print "ERROR - {}: no such file".format(config_file_name)
			exit()
	except IndexError:
		print "\nUSAGE: python {} <quench.cfg>\n".format(sys.argv[0])
		exit()

	os.chdir("../scratch")
	os.system("""
	ln -sf {}/potentials/pseudo.H_pbe
	ln -sf {}/potentials/pseudo.C_pbe
	ln -sf {}/potentials/pseudo.O_pbe
	ln -sf {}/potentials/pseudo.N_soft_pbe
	""".format(maindir, maindir, maindir, maindir))

	for file_name in file_names:
		print "\n----------------------------------------------------------"
		print "INPUT: {}.pdb".format(file_name)
		print "----------------------------------------------------------"
		print "\ngenerating .xyz file for input .pdb"
		command = "python {}/util/pdb2xyz.py {}.pdb > {}.xyz".format(maindir, file_name, file_name)
		print command
		os.system(command)

		print "\nrunning quench"
		command = "srun -ppdebug -n108 {}/bin/mgmol-pel -c ../DNA/{} -i {}.xyz > {}.out".format(maindir, config_file_name, file_name, file_name)
		print command
		os.system(command)

		print "\ncreating .xyz with MLWCs from quench output"
		command = "python {}/util/getMLWCinXYZ.py {}.out > {}-mlwc.xyz".format(maindir, file_name, file_name)
		print command
		os.system(command)

	os.chdir("..")
