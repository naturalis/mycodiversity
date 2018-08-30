import functions
configfile: "config.yml"
###TODO ALL###
# - Add protected and tmp indications


### Rule all creates output files and assigns sample and output directory wildcards ###
rule all:
    input:
        expand("{outdir}/ZOTUS/{sample}_zotutab.txt", sample=config["samples"], outdir=config["outdir"]),
        expand("{outdir}/qual_filter/{sample}_EE_report.txt", sample=config["samples"], outdir=config["outdir"]),
        expand("{outdir}/ZOTUS/abundant/{sample}_zotutab_af.txt", sample=config["samples"], outdir=config["outdir"]),
        expand("{outdir}/FINAL/{sample}_zotus_final.fa", sample=config["samples"], outdir=config["outdir"]),
        expand("{outdir}/FINAL/{sample}_zotutab_final.txt", sample=config["samples"], outdir=config["outdir"])

###  Rule cuts the primers off the reads. If Illumina is used, the rule cuts paired-end. ###
###  \\TODO: minlen should be implemented, as well as unfiltered deletion ###
rule filter_primers:
    input:
        lambda wildcards: functions.getReads(config["samples"][wildcards.sample],config["params"]["platform"],1)
    params:
        FWD=lambda wildcards: config["primers"]["FWD"],
        REV=lambda wildcards: config["primers"]["REV"], # REVCOMP!!
        FWD_RC=functions.getRevComp(config["primers"]["FWD"]),
        REV_RC=functions.getRevComp(config["primers"]["REV"]),
        minlen=functions.changeDefault(config, "cutadapt_minlen","-m 100"),
        platform=lambda wildcards: config["params"]["platform"],
        reversereads=lambda wildcards: functions.getReads(config["samples"][wildcards.sample],config["params"]["platform"],2)
    output:
        R1_out="{outdir}/filtered/{sample}_R1_filt.fastq",
        R2_out="{outdir}/filtered/{sample}_R2_filt.fastq"
    log:
        "{outdir}/log/filter_primers/{sample}_error.log"
    run:
        print({params.platform}, params.platform)
        if params.platform == "illumina":
            shell("cutadapt -g {params.FWD} -a {params.REV_RC} -G {params.REV} -A {params.FWD_RC} -o {output.R1_out} {params.minlen} \
				  -p {output.R2_out} --discard-untrimmed --pair-filter=both {input} {params.reversereads} 2>{log}")
        elif str(params.platform) == "454":
            shell("cutadapt -g {params.FWD} -a {params.REV_RC} -o {output.R1_out} {params.minlen} {input} 2>{log}")
            shell("touch {output.R2_out}")

### If the platform is Illumina, rule calls scripts to decide whether or not to merge the fwd and rev reads ###
### If the platform is not illumina and the amplicon is not the full length ITS region, the reads are truncated ###
### To 250 bp. If amplicon is full length ITS, nothing is done ###
### LEG:         reversereads=lambda wildcards: functions.getReads(config["samples"][wildcards.sample],config["params"]["platform"],2,"filt"),

rule merge_reads:
    input:
        R1="{outdir}/filtered/{sample}_R1_filt.fastq",
    output:
        "{outdir}/merged/{sample}_merged_reads.fq"
    log:
        "{outdir}/log/merged/{sample}_error.log"
    params:
        min_overlap=functions.changeDefault(config, "flash_min_overlap"), # optional param (but might to set a default later //TODO)
        platform=config["params"]["platform"],
        amplicon=config["primers"]["amplicon"],
        outdir=config["outdir"],
        revreads="{outdir}/filtered/{sample}_R2_filt.fastq"
    run:
        if params.platform == "illumina":
            mergebool= functions.decideMerger(wildcards.sample, wildcards.outdir)
            if "MERGE" in mergebool:
                shellstring = "flash {input_r1} {input_r2} -c > {output} {overlap} -M 250 2>{log}".format(input_r1=input.R1, input_r2 = params.revreads, output=output, overlap=params.min_overlap, log=log)
                #shell("flash {input.R1} revreads -c > {output} {params.min_overlap} -M 250 2>{log}")		
                print(shellstring)
                shell(shellstring)
            elif "FWD" in mergebool:
                print("Your reverse reads are bad, I will only use your forward reads\nPlease note that the file name will still be 'merged' (\\\TODO v.03)")
                shell("cp {input.R1} {output}")
                shell("echo Your reverse reads were not used, no merging was done. Please see the documentation for more information > {log}")
        else:
            if params.amplicon != "Full":
                print("Going to truncate")
                shell("./deps/usearch11 -fastx_truncate {input.R1} -trunclen 250 -fastqout {output}")
            else:
                shell("cp {input.R1} {output}")


### Rule creates the USEARCH Estimated Errorr report ###
rule calculate_EE:
    input:
        "{outdir}/merged/{sample}_merged_reads.fq"
    output:
        "{outdir}/qual_filter/{sample}_EE_report.txt"
    log:
        "{outdir}/log/qual_filter/{sample}_EE_error.log"
    shell:
        "./deps/usearch11 -fastx_info {input} -output {output} 2>{log}"

### Rule checks if the user provided an EE threshold, if not it uses the mean value ###
### calculated by the previous rule ###
rule qual_filter:
    input:
        reads="{outdir}/merged/{sample}_merged_reads.fq",
        report="{outdir}/qual_filter/{sample}_EE_report.txt"
    output:
        "{outdir}/qual_filter/{sample}_filtered_reads.fa"
    log:
        "{outdir}/log/qual_filter/{sample}_error.log"
    run:
        max_EE = functions.determineEE(wildcards,config)
        shellstring="./deps/usearch11 -fastq_filter {input.reads} -fastq_maxee %s -fastaout {output} -relabel Filt 2> {log}"%max_EE
        shell(shellstring)


### Rule calls the dereplication function of VSEARCH (not USEARCH due to the 32-bit restrictions) ### 
rule dereplicate:
    input:
        "{outdir}/qual_filter/{sample}_filtered_reads.fa"
    output:
        "{outdir}/dereplicated/{sample}_dereplicated_reads.fa"
    log:
        "{outdir}/log/dereplicate/{sample}_error.log"
    shell:
        "vsearch --derep_fulllength {input} -sizeout -relabel Uniq --output {output} 2> {log}"

### Rule discards singletons ###
rule discard_singletons:
	input:
		"{outdir}/dereplicated/{sample}_dereplicated_reads.fa"
	output:
		"{outdir}/dereplicated/{sample}_desingled_reads.fa"
	log:
		"{outdir}/log/discard_singletons/{sample}_error.log"
	shell:
		"./deps/usearch11 -sortbysize {input} -fastaout {output} -minsize 2 2>{log}"

### Rule runs UNOISE3 algorithm with an alpha of 2 ###
rule denoise:
    input:
        "{outdir}/dereplicated/{sample}_dereplicated_reads.fa"
    output:
        "{outdir}/ZOTUS/{sample}_zotus.fa"
    log:
        "{outdir}/log/denoise/{sample}_error.log"
    params:
        unoise_alpha=functions.changeDefault(config, "unoise_alpha") # optional param
    shell:
        "./deps/usearch11 -unoise3 {input} -zotus {output} -unoise_alpha 2 -minsize 2 2> {log}"

### Rule uses the raw reads to create the ZOTU table indicating the sequential abundance of the ZOTUS ###
rule zotutable:
    input:
        zotus="{outdir}/ZOTUS/{sample}_zotus.fa",
        reads="{outdir}/merged/{sample}_merged_reads.fq"
    output:
        "{outdir}/ZOTUS/{sample}_zotutab.txt"
    log:
        "{outdir}/log/zotutable/{sample}_error.log"
    shell:
        "./deps/usearch11 -otutab {input.reads} -zotus {input.zotus} -otutabout {output} 2>{log}"

### Rule calls the abundance filter script that filters the ZOTU set for 0.5% of the total reads ###
### It then filters out the non-abundant ZOTUS ###
rule abundance_filter:
    input:
        zotus="{outdir}/ZOTUS/{sample}_zotus.fa",
        zotutab="{outdir}/ZOTUS/{sample}_zotutab.txt"
    output:
        af_zotus="{outdir}/ZOTUS/abundant/{sample}_zotus_af.fa",
        af_zotutab="{outdir}/ZOTUS/abundant/{sample}_zotutab_af.txt",
        intermediate="{outdir}/ZOTUS/abundant/{sample}_names.tmp"
    log:
       "{outdir}/log/abundance/{sample}_error.log"
    run:
       shell("python abundance_filter.py {input.zotutab} {output.af_zotutab} {output.intermediate} 0.005 2>{log}")
       shell("./deps/faSomeRecords {input.zotus} {output.intermediate} {output.af_zotus} 2>>{log}") 

### Rule calls the contamination filtering script. This blasts the ZOTUS against the UNITE ###
### database and selects only the ZOTUS that 1) have a hit in UNITE and 2) are matched with###
### a fungus. ###
rule filter_contamination_filter:
    input:
        "{outdir}/ZOTUS/abundant/{sample}_zotus_af.fa"
    output:
        "{outdir}/FINAL/{sample}_zotus_final.fa"
    log:
        '{outdir}/log/contamination/{sample}_error.log'
    shell:
        "bash excludeContamination.sh {input} {wildcards.sample}_blast.txt {output}"


rule filter_zotutab:
    input:
        zotus="{outdir}/FINAL/{sample}_zotus_final.fa",
        zotutab="{outdir}/ZOTUS/abundant/{sample}_zotutab_af.txt"
    output:
        "{outdir}/FINAL/{sample}_zotutab_final.txt"
    log:
        "{outdir}/log/ZOTUS/{sample}_error.log"
    shell:
        """
        grep ">" {input.zotus} | sed 's/>//g' > tmp
        head -1 {input.zotutab} > {output}
        grep -f tmp {input.zotutab} >> {output}
        rm tmp
        """
