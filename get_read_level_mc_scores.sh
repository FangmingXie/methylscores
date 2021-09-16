#!/bin/bash

mydir=$(dirname $0)
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
	echo "Parameter 5 is empty; choose from single, paired, and cgonly"
    exit 1
elif [[ -f "$3" ]]; then
	echo "$3 exists..."
    exit 1
elif [[ -f "$4" ]]; then
	echo "$4 exists..."
    exit 1
fi


if [ "$5" != "single" ] && [ "$5" != "paired" ] && [ "$5" != "cgonly" ]; then
	echo "choose from single, paired, and cgonly"
    exit 1
fi



# 0-based coords
bam=$1
bed=$2
output_mcinfo=$3
output_mcscores=$4
mod=$5

echo $mod
if [[ "$mod" == "single" ]]; then
    echo "collecting info from $bam"
    $mydir/get_mc_pattern_bed.sh $bam $bed ${output_mcinfo} "normal"
    echo "calculating scores from ${output_mcinfo}"
    $mydir/compute_scores.py -i ${output_mcinfo} -o ${output_mcscores}

elif [[ "$mod" == "cgonly" ]]; then
    echo "collecting info from $bam"
    $mydir/get_mc_pattern_bed.sh $bam $bed ${output_mcinfo} "cgonly"
    echo "calculating scores from ${output_mcinfo}"
    $mydir/compute_scores_cgonly.py -i ${output_mcinfo} -o ${output_mcscores} # use cgonly

elif [[ "$mod" = "paired" ]]; then
    echo "collecting info from $bam"
    $mydir/get_mc_pattern_bed.sh $bam $bed ${output_mcinfo} "paired"
    echo "calculating scores from ${output_mcinfo}"
    $mydir/compute_scores.py -i ${output_mcinfo} -o ${output_mcscores}

fi

