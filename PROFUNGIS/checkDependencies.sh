#!/bin/bash

#Colour codes
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'

echo "This script will check if all dependencies for PROFUNGIS are in place"

function checkIfExists {
if [ ! -f ./deps/${1} ]; then
    echo -e "${RED}"
    echo "file ${1} does not exist"
    echo -e "${NC}"
    exit 1
fi
}

toCheck=("flash" "cutadapt" "faSomeRecords" "usearch11" "fasterq-dump"\
 "vsearch" "primer.data" "trimmomatic.jar" "Unite/unite.nhr" "Unite/unite.nin" \
"Unite/unite.nsq")

for tool in ${toCheck[@]}; do
	checkIfExists ${tool}
done


hash blastn 2>/dev/null || { echo -e "${RED}"; echo >&2 "blastn is not installed"; echo -e "${NC}"; exit 1;}
hash snakemake 2>/dev/null || { echo -e "${RED}"; echo >&2 "snakemake is not installed"; echo -e "${NC}"; exit 1;}
hash python3 2>/dev/null || { echo -e "${RED}"; echo >&2 "python3 is not installed"; echo -e "${NC}"; exit 1;}


echo -e "${GREEN}"
echo "All dependencies are in place, you can go on and run the pipeline!"
echo -e "${NC}"
