#!/bin/bash

# set the number of nodes
#SBATCH --nodes=1

# set max wallclock time
#SBATCH --time=12:00:00

# set name of job
#SBATCH --job-name=Bw

# mail alert at start,/data/projects/p283_rna_and_disease/projects/sinergia/samples/fastq/ end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=carlos.pulido@dbmr.unibe.ch

# set memory [K|M|G|T]
#SBATCH --mem=20G

# set working directory "--chdir=./"
#SBATCH

# set error file
#SBATCH --error="%A.error"

# set output file
#SBATCH --output="%A.output"

#extract spliced reads
#file="cFb1Aligned.sortedByCoord.out.bam"
#file="cFb1aAligned.sortedByCoord.out.bam"
#file="cFb2Aligned.sortedByCoord.out.bam"
#file="cFb3Aligned.sortedByCoord.out.bam"
#file="cFb3aAligned.sortedByCoord.out.bam"
#file="iCPC10Aligned.sortedByCoord.out.bam"
#file="iCPC13Aligned.sortedByCoord.out.bam"
#file="iCPC17Aligned.sortedByCoord.out.bam"
#file="iCPC19Aligned.sortedByCoord.out.bam"
file="iCPC8Aligned.sortedByCoord.out.bam"

samtools view -H star/map/$file >> assembly/statistics/$file.sam
samtools view -h star/map/$file | awk '$6 ~ /N/' >> assembly/statistics/$file.sam
samtools view -h -b assembly/statistics/$file.sam > assembly/statistics/$file.spliced.reads.bam
rm assembly/statistics/$file.sam

#index bam file
samtools index -b assembly/statistics/$file.spliced.reads.bam

#convert to bigwig
module add UHTS/Analysis/deepTools/2.5.4
bamCoverage --filterRNAstrand forward --bam assembly/statistics/$file.spliced.reads.bam --outFileFormat bigwig --outFileName assembly/statistics/$file.spliced.reads.forward.bw
bamCoverage --filterRNAstrand reverse --bam assembly/statistics/$file.spliced.reads.bam --outFileFormat bigwig --outFileName assembly/statistics/$file.spliced.reads.reverse.bw
