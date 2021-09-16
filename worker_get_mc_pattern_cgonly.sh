#!/bin/bash
# getting mc patterns for each read overlapping a given query region from bismark bam file
# can choose from short or long output format
# example usage: 
# 01.1.get_mc_pattern.sh test_sample_R1_bismark_bt2_pe.deduplicated.sorted.bam chr11 18279557 18279560

# !!! methylation string column could be $14 or $15

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
elif [[ -z "$4" ]]; then
	echo "Parameter 4 is empty"
    exit 1
fi

if [[ -z "$5" ]]; then
	outputFormat='short' # default short format	
else
	outputFormat=$5 # choose a short or long 
fi

# 1-based coords
bam=$1
chr=$2 # "chr7"
start=$3
end=$4

### prototype1 
## output bed format (0-based coords)
# samtools view $bam "$chr:$start-$end" \
# 	| awk 'BEGIN{FS=OFS="\t"} 
# 			{gsub("XM:Z:", "", $15); 
# 			print $3,$4-1,$4-1+length($15),$15}'

## calculations
## global 0-based coords:
## read ($4-1, $4-1+len)
## query (start-1, end) 

## local 0-based coords:
## read (0, len)
## query (start-$4, end-$4+1)

## to keep: 
## read (lstart, lend) = (max(0, start-$4), min(len, end-$4+1))
## 1-based substr
## substr(s, 
##		lstart+1, 
##		lend-lstart)


if [ "$outputFormat" == "long" ]; then
	### prototype2 
	# chr, start, end, string [0-based, global coords, overlapping with query regions only]
	samtools view -F 4 $bam "$chr:$start-$end" \
		| awk -v start=$start -v end=$end \
			'BEGIN{FS=OFS="\t"} {
				gsub("XM:Z:", "", $15); 
				len=length($15);
				lstart=0;
				if (lstart<start-$4) {lstart=start-$4};
				lend=len;
				if (lend>end-$4+1) {lend=end-$4+1};
				print $3, $4-1+lstart, $4-1+lend, substr($15, lstart+1, lend-lstart)
			}'
elif [ "$outputFormat" == "short" ]; then
	### prototype3
	# chr, querystart, queryend, string, string, string, ...
	samtools view -F 4 $bam "$chr:$start-$end" \
		| awk -v chr=$chr -v start=$start -v end=$end \
			'BEGIN{FS=OFS="\t"; printf "%s\t%d\t%d\t", chr, start-1, end} 
			{
				gsub("XM:Z:", "", $15); 
				len=length($15);
				lstart=0;
				if (lstart<start-$4) {lstart=start-$4};
				lend=len;
				if (lend>end-$4+1) {lend=end-$4+1};
				preout=substr($15, lstart+1, lend-lstart);
				gsub(/\./, "", preout);
				gsub(/h/, "", preout);
				gsub(/H/, "", preout);
				gsub(/x/, "", preout);
				gsub(/X/, "", preout);
				printf "%s,", preout;
			}
			 END{printf "\n"}
			'
else
	echo "choose from short or long!"
fi
