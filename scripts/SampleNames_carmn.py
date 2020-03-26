import sys, os, json

with open("./config/config.json", "r") as read_file:
    config = json.load(read_file)

if not os.path.isdir(config["fastqFolder"]):
	sys.exit("The directory with the fastq files does not exist")
if not os.path.isfile(config["gtf"]):
	sys.exit("File {} does not exist".format(config["gtf"]))
if not os.path.isfile(config["experimental_design"]):
	sys.exit("File {} does not exist".format(config["gtf"]))

SampleIDs = []
SampleNames = []
SamplesPerID = {}

for file in os.listdir(config["fastqFolder"]):
	if file.endswith(config["extension"]):

		ID = file.split(config["delim"])[0]

		name = file.split(".")[0]
		#name = config["delim"].join(name)
		
		SampleIDs.append(ID+config["delim"])
		SampleNames.append(name)

		if not ID in SamplesPerID:
			SamplesPerID[ID] = []

		if not name in SamplesPerID[ID]:
			SamplesPerID[ID].append(name)

SampleIDs = sorted(list(set(SampleIDs)))		
SampleNames = sorted(list(set(SampleNames)))

mates=sorted(list(config["mates"].values()))

fastqSamples = [SampleNames, SampleIDs, SamplesPerID, mates]

for name in SamplesPerID:

	mates1 = []
	mates2 = []
	
	for file in SamplesPerID[name]:
		if "_R1_" in file:
			mates1.append(config["fastqFolder"]+"/"+file+".fastq.gz")
		if "_R2_" in file:
			mates2.append(config["fastqFolder"]+"/"+file+".fastq.gz")

	mates1 = sorted(mates1)
	mates2 = sorted(mates2)

	mates1 = ",".join(mates1)
	mates2 = ",".join(mates2)

	ref = "/data/projects/p283_rna_and_disease/DATA.FILES/STAR_ref/hg38/"
	command = "STAR --outSAMstrandField intronMotif --limitBAMsortRAM 25000000000 --runThreadN 2 --outSAMtype BAM SortedByCoordinate --genomeDir "+ref+" --readFilesCommand zcat --readFilesIn "+mates1+" "+mates2+" --outFileNamePrefix star/map/"+name

	print command