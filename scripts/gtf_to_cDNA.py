# Carlos P. 4/12/2017

# Reads a gtf file and capture specific exons for each transcript
# Use exonic coordinates to extract exon sequences using the providen fasta file
# Prints cDNA sequence for each transcript in fasta format

# USAGE:
# python gtf_to_cDNA.py input.gtf input.fasta > output.fasta 

import sys
from Bio import SeqIO

gtf = sys.argv[1]
fasta = sys.argv[2]
transcripts = {}

fi = open(gtf, 'r')
for line in fi:

	line = line.strip().split("\t")
	flags = line[-1].split(" ")

	if len(line) < 2: continue
	if line[2] == "exon":

		gen_id = flags[1].replace("\"", "").replace(";","")
		transcript_id = flags[3].replace("\"", "").replace(";","")
		chrom = line[0]
		start = line[3]
		stop = line[4]

		if transcript_id in transcripts:

			transcripts[transcript_id].append([chrom, start, stop, gen_id, transcript_id])

		else:

			transcripts[transcript_id] = [[chrom, start, stop, gen_id, transcript_id]]

fi.close()


fasta_sequences = SeqIO.parse(open(fasta),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq)
	
	
	for tr in transcripts:
		tr_seq = ""
		tr_start = ""
		if transcripts[tr][0][0] == name:
			for exon in transcripts[tr]:
				if tr_start == "":
					tr_start = exon[1]
				start = int(exon[1])
				end = int(exon[2])
				tr_seq = tr_seq+sequence[start-1:end]
			print ">"+exon[3]+";"+tr+"|"+name+":"+tr_start+"-"+exon[2]+"\n"+tr_seq

