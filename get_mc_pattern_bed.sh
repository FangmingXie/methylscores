#!/bin/bash

# given multiple regions (bed file), getting mc patterns for each read overlapping a given query region from bismark bam file
# This is a wrap up of the 01.1.get_mc_pattern.sh
mydir=$(dirname $0)

# check parameters
if [[ -z "$1" ]]; then
    echo "Parameter 1 is empty"
    exit 1
elif [[ -z "$2" ]]; then
	echo "Parameter 2 is empty"
    exit 1
elif [[ -z "$3" ]]; then
	echo "Parameter 3 is empty"
    exit 1
elif [[ -f "$3" ]]; then
	echo "$3 exists..."
    exit 1
fi


# 0-based coords
bam=$1
bed=$2
output=$3

# read through bed file
while read p || [ -n "$p" ]
do
	row=(${p//"\t"/" "})
	chr=${row[0]}
	start=${row[1]}
	end=${row[2]}
	$mydir/worker_get_mc_pattern.sh $bam $chr $((start+1)) $end >> $output
done < $bed

