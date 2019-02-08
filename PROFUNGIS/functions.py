import subprocess
import sys

def numSeqs(file):
	"""\
	Function determines the number of sequences in the input file using egrep
	"""
	file = str(file[0])
	proc=subprocess.Popen("egrep -c '>' %s"%(file), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	numseqs = proc.communicate()
	numseqs = numseqs[0].decode("utf-8").strip("\n")

	return(int(numseqs))

def getReads(sample, platform,orient):
	"""\
	Function determines whether or not the sample has a reverse read file.
	If the platform is illumina, the reads are returned. If the platform
	is a pyrosequencer, there is no reverse read file, and the only file is returned.
	"""
	pyros = ["454","iontorrent"]
	if str(platform) in pyros and orient == 1:
		return("{sample}.fastq".format(sample=sample))
	if platform == "illumina":
		return("{sample}_{orient}.fastq".format(sample=sample, orient=orient))
	else:
		return("")

def getRevComp(primer):
	"""\
	Function determines the reverse complement of a string. The complement
	dictionary contains all possibilities of bases/wildcards and their
	complement.
	"""
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'R':'Y','Y':'R',
				  'S':'S', 'W':'W', 'K':'M','B':'T','V':'B', 'D':'H','H':'D'}
	reverse_complement = "".join(
		complement.get(base, base) for base in reversed(primer)
		)
	return(reverse_complement)

def changeDefault(config, param, default=""):
    """\
    Checks if a default value is given. If not, the parameter is optional
    and both the flag and the value need to be set. If the paramter is
    not set in the config file, the parameter is not used
    If default is set, the parameter is obligatory and if there is no
    user set value in the config file, the default is used.
    """
    try:
        if default == "":
            outstring="%s %s"%(config["params"][param]["flag"],
			config["params"][param]["value"])
        else:
            outstring="%s %s"%(config["params"][param]["flag"],
			config["params"][param]["value"])
        return(outstring)
    except KeyError:
        return(default)


def determineEE(wildcards, config):
	"""\
	Function determines the calculated Estimated Error (EE) scores
	calculated by UNOISE. It first checks if there is a value entered
	in the config file. If not, the output of UNOISE is parsed
	to a dictionary containing all EE stats. By default,
	the function returns the mean value but, this can be
	altered to be the other values provided by UNOISE.
	If the config file contains an entry for max_EE,
	this value is chosen and the user is informed so.
		This function is currently unused as the default value is automatically
		used (EE = 1.0)
	"""
	STAT='mean' # CHANGE ME
	value= changeDefault(config, "usearch_max_EE", "not set")
	sample = wildcards.sample
	outdir = wildcards.outdir
	if value == "not set":
		report=open("./%s/qual_filter/%s_EE_report.txt"%(outdir,sample), 'r').readlines()
		EEstats = report[-1].replace(";",",").replace("EE","")
		statlist = EEstats.split(',')
		statlist = [x.strip() for x in statlist]
		statdict = dict([x.split(' ') for x in statlist])
		return(statdict[STAT])
	else:
		print("Using user input for max EE threshold of %s for sample %s"%(value, sample))
		return(value)

def decideMerger(sample, outdir):
    """\
    Function evaluates the decideMerger script. The forward and reverse samples
	are determined, and the script is called with these read files.
	The output of that script is parsed to subprocess.PIPE and this is parsed
	into a string that is returned.
    """
    fwd="%s/filtered/%s_R1_filt.fastq"%(outdir, sample)
    rev="%s/filtered/%s_R2_filt.fastq"%(outdir, sample)
    proc=subprocess.Popen("bash decideMerger.sh %s %s"%(fwd,rev) , shell=True,
						   stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    decision = proc.communicate()
    decision = decision[0].decode("utf-8").strip("\n")
    print(decision)
    return(decision)
