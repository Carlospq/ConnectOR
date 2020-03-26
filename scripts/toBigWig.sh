#!/bin/bash

#$1	deeptools version
#$2	samtools versoin
#$3	sample name
#$4	type of library (single/paired)
#$5	project name

deeptools=$1
samtools=$2
sample_name=$3
project_name=$4
n=$5

module add UHTS/Analysis/deepTools/$deeptools
module add UHTS/Analysis/samtools/$samtools
mkdir -p $project_name/star/bigwig

extension="Aligned.sortedByCoord.out.bam"

#Generate index from BAM
samtools index -b $project_name/star/map/$sample_name$extension

#FORWARD reads
bamCoverage -b $project_name/star/map/$sample_name$extension \
 -o $project_name/star/bigwig/$sample_name.fwd.bw \
 --filterRNAstrand forward \
 --normalizeUsingRPKM \
 --samFlagInclude 64 \
 --numberOfProcessors $n

#REVERSE reads
bamCoverage -b $project_name/star/map/$sample_name$extension \
 -o $project_name/star/bigwig/$sample_name.rev.bw \
 --filterRNAstrand reverse \
 --normalizeUsingRPKM \
 --samFlagInclude 64 \
 --numberOfProcessors $n
