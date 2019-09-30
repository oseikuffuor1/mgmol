import os
import re
import sys

# This script runs quench on each tripeptide by doing slight domain alterations to the provided .cfg file
# run it from the bio-mlwc/protein directory
# make sure the bio-mlwc/scratch directory exists so the script can move there to run the computations

if __name__ == "__main__":
	# neutral residues
	# make sure the charge is 0 in the .cfg file
	residues = [
	"ALA", "ASN", "CYS", "GLN",
	"GLY", "VAL", "ILE", "LEU",
	"MET", "PHE", "PRO", "SER",
	"THR", "TRP", "TYR", "HSP",
	"HSD", "HSE"
	]

	# # acidic residues
	# # make sure the charge is -1 in the .cfg file
	# residues = [
	# "ASP", "GLU"
	# ]

	# # basic residues
	# # make sure the charge is 1 in the .cfg file
	# residues = [
	# "ARG", "LYS"
	# ]

	file_names = ["../protein/template_files/{}/{}-tripeptide".format(r, r) for r in residues]

	maindir = "../../.."

	try:
		config_file_name = sys.argv[1]
        	with open(config_file_name, "r") as config_file:
                	config = config_file.read()
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

	bufsize=6
	ox_pattern = re.compile("ox=-*[0-9]+\.")
	oy_pattern = re.compile("oy=-*[0-9]+\.")
	oz_pattern = re.compile("oz=-*[0-9]+\.")

	for file_name in file_names:
		with open("{}.pdb".format(file_name), "r") as f:
			lines = f.readlines()

		ang2bohr=1.8897269

		minx=1000
		miny=1000
		minz=1000
		for line in lines:
		  if 'ATOM' in line[0:6] or 'HETATM' in line[0:6]:
			x = float(line[30:38])
			y = float(line[38:46])
			z = float(line[46:54])
			minx=min(x,minx)
			miny=min(y,miny)
			minz=min(z,minz)

		ox, oy, oz = ang2bohr*minx-bufsize, ang2bohr*miny-bufsize, ang2bohr*minz-bufsize

		new_config = config
		new_config = re.sub(ox_pattern, "ox={}.".format(str(int(round(ox)))), new_config)
		new_config = re.sub(oy_pattern, "oy={}.".format(str(int(round(oy)))), new_config)
		new_config = re.sub(oz_pattern, "oz={}.".format(str(int(round(oz)))), new_config)

		with open("{}.cfg".format(file_name), "w") as new_config_file:
			new_config_file.write(new_config)

		print "\n----------------------------------------------------------"
		print "INPUT: {}.pdb".format(file_name)
		print "----------------------------------------------------------"
		print "\ngenerating .xyz file for input .pdb"
		command = "python {}/util/pdb2xyz.py {}.pdb > {}.xyz".format(maindir, file_name, file_name)
		print command
		os.system(command)

		print "\nrunning quench"
		command = "srun -ppdebug -n48 {}/bin/mgmol-pel -c {}.cfg -i {}.xyz > {}.out".format(maindir, file_name, file_name, file_name)
		print command
		os.system(command)

		print "\ncreating .xyz with MLWCs from quench output"
		command = "python {}/util/getMLWCinXYZ.py {}.out > {}-mlwc.xyz".format(maindir, file_name, file_name)
		print command
		os.system(command)

	os.chdir("..")
