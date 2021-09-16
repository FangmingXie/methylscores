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
## to keep global
## read (gstart, gend) = (max($4-1, start-1), min($4-1+len, end))

## local 0-based coords (-($4-1)):
## read (0, len)
## query (start-$4, end-$4+1)
## to keep local 
## read (lstart, lend) = (max(0, start-$4), min(len, end-$4+1))
## 1-based substr: substr(s, lstart+1, lend-lstart)

## 99 (read + strand; R1); 147 (read - strand; R2) mate
## 83 (read - strand; R1); 163 (read + strand; R2) mate

# need to align reads in pairs
if [ "$outputFormat" == "test" ]; then
	### prototype2 
	# chr, start, end, string [0-based, global coords, overlapping with query regions only]
	samtools view -F 4 -f 2 $bam "$chr:$start-$end" \
		| sort -k 1,1 \
		| awk -v start=$start -v end=$end \
		'BEGIN{FS=OFS="\t" 
			   lastReadName=""
			   lastSeq="."
			  } 
			{
				# readname $1, chr $3, start $4, sequence $15
				gsub("XM:Z:", "", $15) 
				len=length($15)
				# max
				lstart=0
				if (lstart<start-$4) {lstart=start-$4}
				# min
				lend=len
				if (lend>end-$4+1) {lend=end-$4+1}

				# test if read name is the same
				if (lastReadName==$1) {
					# paired with last row
					# should combine this row with the last row
					# should not release this row, it will be print next time

					currStart=$4-1+lstart
					currEnd=$4-1+lend
					currSeq=substr($15, lstart+1, lend-lstart)

					# get min max; smallEnd, largeStart
					if (lastStart < currStart) {
						smallStart = lastStart
						smallEnd = lastEnd
						smallSeq = lastSeq
						largeStart = currStart
						largeEnd = currEnd
						largeSeq = currSeq
					}
					else {
						smallStart = currStart
						smallEnd = currEnd
						smallSeq = currSeq
						largeStart = lastStart
						largeEnd = lastEnd
						largeSeq = lastSeq
					}
					# truncate
					trStart = 0
					trEnd = length(largeSeq) 
					if (smallEnd > largeStart) {
						trStart = smallEnd - largeStart
					}

					combinedSeq = (smallSeq ";" substr(largeSeq, trStart+1, trEnd-trStart))

					# update 
					lastStart=(smallStart ";" largeStart) # does not matter
					lastEnd=(smallEnd ";" largeEnd) # does not matter
					lastSeq=combinedSeq
				}
				else {
					# not paired with last row
					# should release last row
					print lastChr, lastStart, lastEnd, lastSeq

					# update 
					lastReadName=$1
					lastChr=$3
					lastStart=$4-1+lstart
					lastEnd=$4-1+lend
					lastSeq=substr($15, lstart+1, lend-lstart)
				}

			}
		END{print lastChr, lastStart, lastEnd, lastSeq}	
		'
elif [ "$outputFormat" == "long" ]; then
	### prototype2 not maintained 
	# chr, start, end, string [0-based, global coords, overlapping with query regions only]
	samtools view -F 4 -f 2 $bam "$chr:$start-$end" \
		| sort -k 1,1 \
		| awk -v start=$start -v end=$end \
		'BEGIN{FS=OFS="\t" 
			   lastReadName=""
			   lastSeq="."
			  } 
			{
				# readname $1, chr $3, start $4, sequence $15
				
				gsub("XM:Z:", "", $15) 
				len=length($15)
				# max
				lstart=0
				if (lstart<start-$4) {lstart=start-$4}
				# min
				lend=len
				if (lend>end-$4+1) {lend=end-$4+1}

				# test if read name is the same
				if (lastReadName==$1) {
					print "!"$3, $4-1+lstart, $4-1+lend, $1, substr($15, lstart+1, lend-lstart)
				}
				else {
					print $3, $4-1+lstart, $4-1+lend, $1, substr($15, lstart+1, lend-lstart)
				}
				# update 
				lastReadName=$1
			}
		'
elif [ "$outputFormat" == "short" ]; then
	### prototype3
	# chr, querystart, queryend, string, string, string, ...
	samtools view -F 4 -f 2 $bam "$chr:$start-$end" \
		| sort -k 1,1 \
		| awk -v chr=$chr -v start=$start -v end=$end \
		'BEGIN{FS=OFS="\t" 
			   lastReadName=""
			   lastSeq="."
			   printf "%s\t%d\t%d\t", chr, start-1, end
			  } 
			{
				# readname $1, chr $3, start $4, sequence $15
				gsub("XM:Z:", "", $15) 
				len=length($15)
				# max
				lstart=0
				if (lstart<start-$4) {lstart=start-$4}
				# min
				lend=len
				if (lend>end-$4+1) {lend=end-$4+1}

				# test if read name is the same
				if (lastReadName==$1) {
					# paired with last row
					# should combine this row with the last row
					# should not release this row, it will be print next time

					currStart=$4-1+lstart
					currEnd=$4-1+lend
					currSeq=substr($15, lstart+1, lend-lstart)

					# get min max; smallEnd, largeStart
					if (lastStart < currStart) {
						smallStart = lastStart
						smallEnd = lastEnd
						smallSeq = lastSeq
						largeStart = currStart
						largeEnd = currEnd
						largeSeq = currSeq
					}
					else {
						smallStart = currStart
						smallEnd = currEnd
						smallSeq = currSeq
						largeStart = lastStart
						largeEnd = lastEnd
						largeSeq = lastSeq
					}
					# truncate
					trStart = 0
					trEnd = length(largeSeq) 
					if (smallEnd > largeStart) {
						trStart = smallEnd - largeStart
					}

					combinedSeq = (smallSeq ";" substr(largeSeq, trStart+1, trEnd-trStart))

					# update 
					lastStart=(smallStart ";" largeStart) # does not matter
					lastEnd=(smallEnd ";" largeEnd) # does not matter
					lastSeq=combinedSeq
				}
				else {
					# not paired with last row
					# should release last row
					printf "%s,", lastSeq 

					# update 
					lastReadName=$1
					lastChr=$3
					lastStart=$4-1+lstart
					lastEnd=$4-1+lend
					lastSeq=substr($15, lstart+1, lend-lstart)
				}

			}
		END{printf "%s,", lastSeq; printf "\n"}
		'
else
	echo "choose from short or long!"
fi
