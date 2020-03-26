#!/bin/bash
if [ ! -s $2 ]; then
	echo Sample\ ID,Number\ of\ input\ reads,Average\ input\ read\ length,Average\ mapped\ length,Mismatch\ rate\ per\ base,Uniquely\ mapped\ reads\ %,Number\ of\ reads\ mapped\ to\ multiple\ loci,%\ of\ reads\ unmapped | awk 'BEGIN {FS = ","}; {OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > $2
fi

flags=( Number\ of\ input\ reads Average\ input\ read\ length Average\ mapped\ length Mismatch\ rate\ per\ base Uniquely\ mapped\ reads\ % Number\ of\ reads\ mapped\ to\ multiple\ loci %\ of\ reads\ unmapped )
for ((i = 0; i < ${#flags[@]}; i++)); do
	if [[ ${flags[$i]} == *'unmapped'* ]]; then
		mm=$( grep %\ of\ reads\ unmapped $1 | awk '{split($NF,a,"%"); sum += a[1]} END {print sum}' )
	else
		IFS="|"
		value=$( grep -E "${flags[$i]}" $1 | awk '{split($NF,a,"%"); print a[1]}' )
		new_line+=$value","
	fi
done
echo $1","$new_line$mm | awk 'BEGIN {FS = ","}; {OFS="\t"; split($1,a,"Log"); split(a[1],b,"/"); print b[length(b)],$2,$3,$4,$5,$6,$7,$8}' >> $2
