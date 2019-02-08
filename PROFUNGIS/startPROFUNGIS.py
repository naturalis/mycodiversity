import argparse
import subprocess
import os.path
import sys
import textwrap
import datetime

from primers import Primer


def normalize(platform):
    """\
    Function takes the platform specificication by the user
    And  checks for synonyms. (Needs updating along the way).
    It returns the platform it recognizes, or quits the program
    when no suitable platform is found.
    """
    platform = platform.lower()
    platforms = ["454","illumina","iontorrent"]
    if platform in platforms:
        return platform
    else:
        print("Platform was not recognized",
              "Please provide either 454, iontorrent or illumina")
        sys.exit()

def paramparser():
    """\
    Function parses commandline arguments and returns all arguments
    and a list with the present optional parameters.
    """
    parser = argparse.ArgumentParser()
    # An SRA or a file with SRA ID's is required, but only one
    id_group = parser.add_mutually_exclusive_group(required=True)
    id_group.add_argument("-m","--multirun", help="Provide file with SRA ids")
    id_group.add_argument("-s","--SRA", type=str, help="Single SRA id")

    parser.add_argument("-f","--forward", required=True, help="Name of forward primer")
    parser.add_argument("-r","--reverse", required=True, help="Name of reverse primer")
    parser.add_argument("-p","--platform",required=True, help="Sequencing platform [illumina|454|iontorrent]")
    parser.add_argument("-l","--local", action="store_true", help="disable SRA downloading, run local")
    # Max Estimated Error for the USEARCH quality filte
    parser.add_argument("-E","--maxEE",help="Estimated Error (EE) filter threshold")
    # Minimal overlap (bp) when merging with FLASH
    parser.add_argument("-o","--minOverlap", help="Minimal overlap when merging reads, 60 by default")
    # Minimal length of the fragment after primer trimming
    parser.add_argument("-L","--minLen", help="Minimal sequence length")
    parser.add_argument("-O","--outdir", help="output directory")
    args = parser.parse_args()

    paramlist = []
    if args.maxEE: paramlist.append(['usearch_max_EE',False,args.maxEE])
    if args.minOverlap: paramlist.append(['flash_min_overlap','-m',args.minOverlap])
    if args.minLen: paramlist.append(['cutadapt_minlen','-m', args.minLen])
    paramlist.append(['platform',False,normalize(args.platform)])
    return args, paramlist

def check_existence(sra_idlist, platform):
    """\
    Function checks if all read files are present in order to run locally
    If the files (or one of the files for illumina) are missing, the pipeline
    will terminate
    """
    errorcount=0
    for id in sra_idlist:
        if platform == "illumina":
            if not os.path.isfile("samples/{sra}/{sra}_1.fastq".format(sra=id)) \
            or not os.path.isfile("samples/{sra}/{sra}_2.fastq".format(sra=id)):
                print("File samples/{sra}/{sra}_1.fastq or samples/{sra}/{sra}_2.fastq does not exist.".format(sra=id))
                errorcount+=1
        elif platform == "454":
            if not os.path.isfile("samples/{sra}/{sra}.fastq".format(sra=id)):
                print("File samples/{sra}/{sra}.fastq does not exist".format(sra=id))
                errorcount +=1
    if errorcount > 0:
        return False
    else:
        return True

def check_amplicon(frag1, frag2):
    """\
    Function checks, based on the primers, which part is sequenced.
    If the target region of the fwd and rev primers are the same, that
    region is used. If they are not the same, the target region is full-length
    ITS. Unless the reverse primer is for ITS1 and the fwd primer is for ITS2,
    which wouldn't be possible, so an error is raised and the pipeline is terminated
    """
    if frag1 == frag2:
        amplicon = frag1
    elif frag1 == "ITS1" and frag2 == "ITS2":
        amplicon = "Full"
    else:
        print("Your forward primer is behind the reverse primer, please check the fragment assignment")
        sys.exit()
    return amplicon

def obtain_primer(primername, primerobject):
    """\
    The script tries to get the sequence attached with the given primer name
    If this name is not present in the primer dataset, the user gets asked
    to add it to the dataset itself via input.
    Next, the user is asked which part of ITS this primer amplifies.
    The primer sequence is returned.
    """
    primersequence = primerobject.get_primer(primername)

    if not primersequence:
        print("Primer {primer} has not yet been used. Please add it if you want\
        to use it".format(primer=primername))
        while True:
            new = input("Input the primer sequence and press enter: ")
            new_frag= input("Input the target sequence of this primer and press\
             enter (ITS1, ITS2): ")
            proper, message = primerobject.add_primer(primername, new, new_frag)
            if not proper:
                print(message)
            else:
                primersequence = new
                fragment = new_frag
                break
    else:
        fragment = primerobject.get_fragment(primername)
    return primersequence, fragment

def create_config(fwd, rev, sra_id, paramlist, outdir ,amplicon):
    """\
    Function creates a yml config file for snakemake. It makes a
    section with all samples, a section with primers and an optional
    section for optional parameters.
    """
    samplestring="samples:\n"
    for id in sra_id:
        samplestring+="  {sample}: samples/{sample}/{sample}\n".format(sample=id)

    outdirstring="outdir: {outdir}".format(outdir=outdir)

    primerstring=textwrap.dedent("""\
    primers:
      FWD: {fwd}
      REV: {rev}
      amplicon: {amplicon}""".format(fwd=fwd, rev=rev, amplicon=amplicon))
    paramstring=""
    if len(paramlist) > 0:
        paramstring=textwrap.dedent("""\
        params:
        """)
        for param in paramlist:
            if not param[1]:
                paramstring+="  {name}: {value}\n".format(name=param[0],value=param[2])
            else:
                paramstring+="  "
                paramstring+=textwrap.dedent("""\
                {name}:
                    flag: {flag}
                    value: {value}
                """.format(name=param[0],flag=param[1],value=param[2]))
    returnstring="{}\n{}\n{}\n{}".format(samplestring,primerstring,paramstring, outdirstring)
    return returnstring

def main():
    """\
    Main function. The arguments and parameters are obtained from
    the parmparser function. A primer object is created and the
    user given names of the primers are attached to the correct sequences.
    The list of SRA ids is determined (>= 1) . Using the SRA-toolkit,
    these ID's are downloaded and stored in the samples folder.
    Next, the config file is created using the samples and the parameters.
    Finally, snakemake is called to continue the pipeline.
    """

    args, paramlist = paramparser()
    primers = Primer()
    fwd,fwd_frag = obtain_primer(args.forward, primers)
    rev, rev_frag = obtain_primer(args.reverse, primers)
    amplicon = check_amplicon(fwd_frag, rev_frag)

    if args.multirun:
        sralist = open(args.multirun,'r').readlines()
        sralist = [x.strip('\n') for x in sralist]
    elif args.SRA:
        sralist = [args.SRA]

    if not args.local:
        print("Start downloading SRA entry")
        for id in sralist:
            os.system("./deps/fasterq-dump {sra} -O ./samples/{sra} --skip-technical".format(sra=id))
            #os.system("./deps/sratoolkit/bin/fasterq-dump {sra} -O ./samples/{sra}".format(sra=id))

    if not check_existence(sralist, args.platform):
        print("One or more of the given SRA entries don't exist locally. Remove the --local flag to download it")
        sys.exit()

    # Next two params need to know what day and time it is
    now = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M")
    if args.outdir:
        outputdir = args.outdir
    else:
        outputdir = "{}_output".format(now)

    config = create_config(fwd,rev,sralist,paramlist, outputdir, amplicon)
    if os.path.isfile("config.yml"):
        timestamp = "config_{}.yml".format(now)
        os.system("cp config.yml {}".format(timestamp))
    with open('config.yml','w') as outfile:
        outfile.write(config)

    subprocess.call(["snakemake"])


if __name__ == "__main__":
    main()
