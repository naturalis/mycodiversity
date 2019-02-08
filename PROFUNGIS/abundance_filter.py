#!/usr/bin/python3
import sys
def main(args):
	"""\
	Function selects the names of ZOTUS with more than x percent of total reads
	assigned. It writes a new OTUtable and a list of ZOTU names
	"""
	otutab_in = args[1]
	otutab_out = args[2]
	otunames_out = args[3]
	otulist = []
	countlist = []
	min_abundance = args[4]
	header = []
	skip = True
	with open(otutab_in,'r') as otutab:
		for line in otutab:
			""" The first line of the OTUtable contains header info.
			It's not needed for the calculation, but used in the
			output file """
			if skip:
				skip = False
				header = line
				continue
			sepline = line.split("\t")
			otulist.append(sepline[0])
			countlist.append(int(sepline[1].strip("\n")))

	# Total number of reads is determined
	total_count = sum(countlist)
	# The threshold is determined based on the cutoff argument
	threshold = total_count*float(min_abundance)

	present_otus = []
	for x in range(len(countlist)):
		if countlist[x] > threshold:
			present_otus.append([otulist[x], countlist[x]])
		else:
			countlist[x] = 0

	outtab = open(otutab_out, 'w')
	outtab.write(header)
	outfa = open(otunames_out,'w')
	for otuset in present_otus:
			outline = str("\t".join(str(x) for x in otuset))
			outtab.write("{}\n".format(outline))
			outfa.write("{}\n".format(otuset[0]))
	outfa.close()
	outtab.close()

if __name__ == "__main__":
	args = sys.argv
	main(args)
