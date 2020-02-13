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

		name = file.split(config["delim"])
		name = config["delim"].join(name[:-1])

		SampleIDs.append(ID)
		SampleNames.append(name+config["delim"])

		if not ID in SamplesPerID:
			SamplesPerID[ID] = []

		if not name+config["delim"] in SamplesPerID[ID]:
			SamplesPerID[ID].append(name+config["delim"])

SampleIDs = sorted(list(set(SampleIDs)))		
SampleNames = sorted(list(set(SampleNames)))

new_SamplesPerID = {}
for ID in SampleIDs:
	new_SamplesPerID[ID] = SamplesPerID[ID]
SamplesPerID = new_SamplesPerID

mates=sorted(list(config["mates"].values()))

fastqSamples = [SampleNames, SampleIDs, SamplesPerID, mates]

groups = {}
for line in open(config["experimental_design"], 'r'):
	line = line.strip().split("\t")
	s = line[0]
	g = line[1]
	if s =="sample": continue
	if not g in groups:
		groups[g] = [s]
	else:
		groups[g].append(s)



