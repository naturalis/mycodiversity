#!/bin/bash

R1=${1}
basename_R1=$(basename ${R1} .fastq)
basename_R2=$(basename ${R2} .fastq)
sample_name=$(basename ${R1} _R1_filt.fastq)
R2=${2}
total_reads=$(egrep -c "@S" ${R1})

### Step 1: do trimmomatic to see if the R2 is worse than R1 ###
java -jar ./deps/trimmomatic.jar SE "${R1}" "${sample_name}_trimmed_1.fastq" SLIDINGWINDOW:4:20 MINLEN:100 2>>log.txt
unfiltered_R1=$(egrep -c "@S" ${sample_name}_trimmed_1.fastq) 
java -jar ./deps/trimmomatic.jar SE "${R2}" "${sample_name}_trimmed_2.fastq" SLIDINGWINDOW:4:20 MINLEN:100 2>>log.txt
unfiltered_R2=$(egrep -c "@S" ${sample_name}_trimmed_2.fastq)

### Step 2 : merge the untrimmed fastq files to see how much can be merged without trimming ###
flash ${R1} ${R2} -m 30 -M 300 -c > ${sample_name}_merged.fastq 2>> log.txt
merged_reads=$(egrep -c "@S" ${sample_name}_merged.fastq)

### Step 3: do trimmomatic on the merged Fastq to see how much gets discarded when merged ###
java -jar ./deps/Trimmomatic/trimmomatic.jar SE "${sample_name}_merged.fastq" "${sample_name}_merged_trimmed.fastq" SLIDINGWINDOW:4:20 MINLEN:10 2>>log.txt
unfiltered_merged=$(egrep -c "@S" ${sample_name}_merged_trimmed.fastq)
#echo "after merging and trimming, ${unfiltered_merged} reads remained"


### Step 4: make a decision ###
### Flags for decisionmaking ###
pairdiff=0
filterflag=0
mergeflag=0
mergefilterflag=0
mergefilterflag2=0

## Step 4.1: see if R1 has better reads than R2 (with a margin of 10%)##
if [ $((${unfiltered_R1}*9/10)) -gt ${unfiltered_R2} ]; then
pairdiff=1
fi

## Step 4.2: see how much of the R2 gets discarded relative to total reads  ##
## If more than 25% of reads is discarded, it gets flagged ##
if [ $((${total_reads}*75/100)) -gt ${unfiltered_R2} ]; then
filterflag=1
fi

## Step 4.3: see how much remain after merging ##
## If more than 25% of reads can't be merged, it gets flagged ##
if [ $((${total_reads}*75/100)) -gt ${merged_reads} ]; then
mergeflag=1
fi

## Step 4.4: see how much remains after filtering the merged reads ##
## Twofold; if more than 25% is discarded, it gets flagged (again) ##
## If the merged-filtered is closer to the original R2 than to the ##
## original R1 percentages, it gets flagged (as the weight of the  ##
## R2 pulls the percentage of good reads down more than the R1     ##
## pulls it up 													   ##

# Step 4.4.1 #
if [ $((${total_reads}/100*75)) -gt ${unfiltered_merged} ]; then
mergefilterflag=1
fi

# Step 4.4.2 uses the formula 2a=< x+y where a is the perc. of reads #
# after merging and x and y are the R1 and R2 percentages. In other  #
# words: is it less than the average of R1 and R2					 #
a=$((${unfiltered_merged}*2))
#echo $((${a})), $((${unfiltered_R1}+ ${unfiltered_R2}))
if [ ${a} -lt $((${unfiltered_R1} + ${unfiltered_R2})) ]; then
mergefilterflag2=1
fi
#echo "${pairdiff}, ${filterflag}, ${mergeflag}, ${mergefilterflag}, ${mergefilterflag2}"
decision=$((${pairdiff}+${filterflag}+${mergeflag}+${mergefilterflag}+${mergefilterflag2}))
if [ ${decision} -gt 2 ]; then
echo "FWD"
else echo "MERGE"
fi


rm "${sample_name}_trimmed_1.fastq"
rm "${sample_name}_trimmed_2.fastq"
rm "${sample_name}_merged.fastq"
rm "${sample_name}_merged_trimmed.fastq"
