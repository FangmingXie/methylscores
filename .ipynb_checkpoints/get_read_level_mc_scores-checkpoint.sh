#!/bin/bash

# check parameters
# bam bed otuput_mc; output_mcscores; mod
if [[ -z "$1" ]]; then
    echo "Parameter 1 is empty"
    exit 1
elif [[ -z "$2" ]]; then
	echo "Parameter 2 is empty"
    exit 1
elif [[ -z "$3" ]]; then
	echo "Parameter 3 is empty"
    exit 1
elif [[ -z "$4" ]]; then
	echo "Parameter 4 is empty"
    exit 1
elif [[ -z "$5" ]]; then
	echo "Parameter 5 is empty; choose from single and paired"
    exit 1
elif [[ -f "$3" ]]; then
	echo "$3 exists..."
    exit 1
elif [[ -f "$4" ]]; then
	echo "$4 exists..."
    exit 1
fi


if [ "$5" != "single" ] && [ "$5" != "paired" ]; then
	echo "choose from single and paired"
    exit 1
fi



# 0-based coords
bam=$1
bed=$2
output_mcinfo=$3
output_mcscores=$4
mod=$5

if [ "$mod" = "single" ]; then
    echo "single-end mode"
    echo "collecting info from $bam"
    01.get_mc_pattern_bed.sh $bam $bed ${output_mcinfo}
    echo "calculating scores from ${output_mcinfo}"
    02.compute_scores.py -i ${output_mcinfo} -o ${output_mcscores}
elif [ "$mod" = "paired" ]; then
    echo "paired-end mode"
    echo "collecting info from $bam"
    01.get_mc_pattern_bed_pairedend.sh $bam $bed ${output_mcinfo}
    echo "calculating scores from ${output_mcinfo}"
    02.compute_scores.py -i ${output_mcinfo} -o ${output_mcscores}
fi

