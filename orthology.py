import os, sys, subprocess
import itertools, gzip

minMatch=0.95

#############################################################################
def read_config(config):
	config_file = open(config, 'r')
	species={}
	for line in config_file:
		#skip header
		if line.startswith("species_name"): continue	
		#generate dictionary with the info for each specie
		line=line.strip().split("\t")
		n,v,p,f=line[0],line[1],line[2],line[3]
		species[n]={"path": p, "version": v, "format": f}
	return(species)

def create_dir(dirName):
	try:
		# Create target Directory
		os.mkdir(dirName)
	except FileExistsError:
		pass

def arguments_dic(line):
	args = line[-1].split(" ")
	dic, n = {}, 0
	for e in args:
		if n == 0:
			key = e
			n = 1
			continue
		if n == 1:
			dic[key] = e.replace("\"","").replace(";","")
			n = 0
	return dic

def retrieve_chainmap_FilesNames(sp1v, sp2v):
	chainmap1=sp1v.lower()+"To"+sp2v.capitalize()+".over.chain.gz"
	chainmap2=sp2v.lower()+"To"+sp1v.capitalize()+".over.chain.gz"
	return([chainmap1,chainmap2])

def retrieve_chainmap_URLs(sp1v, sp2v, chainmaps):
	url_chainmap1='http://hgdownload.cse.ucsc.edu/goldenPath/{}/liftOver/{}'.format(sp1v.lower(), chainmaps[0])
	url_chainmap2='http://hgdownload.cse.ucsc.edu/goldenPath/{}/liftOver/{}'.format(sp2v.lower(), chainmaps[1])
	return([url_chainmap1,url_chainmap2])

def download_chainmaps(urls):
	for url in urls:
		chain=url.split("/")[-1]
		print("Downloading... "+chain, end='\r')
		if not os.path.isfile('./chainmaps/'+chain): subprocess.call(["wget", "-Pchainmaps", url], stderr=open("err.log", 'wb'), shell=False)
		print(' '*14+"\rComplete")

def retrieve_liftover_arguments(species, sp1, sp2, sp1v, sp2v, chainmaps, minMatch):
	arguments=["./liftOver"]
	arguments.append("-minMatch={}".format(minMatch))
	arguments.append(species[sp1]["path"])
	arguments.append("./chainmaps/{}".format(chainmaps[0]))
	arguments.append("{}To{}.{}".format(sp1v,sp2v,"lifover"))
	arguments.append("{}To{}.{}".format(sp1v,sp2v,"unmapped"))
	arguments.append('' if species[sp1]["format"].lower() == "bed" else "-gff")
	return(arguments)

def get_id(line):
	arguments = arguments_dic(line)
	feature = line[2]
	if feature == "gene":
		id_ = "gene:"+arguments["gene_id"]
	elif feature == "transcript":
		id_ = "transcript:"+arguments["gene_id"]+";"+arguments["transcript_id"]
	elif feature == "exon":
		id_ = "exon:"+arguments["gene_id"]+";"+arguments["transcript_id"]+";"+arguments["exon_number"]
	else:
		id_ = feature+":"+arguments["gene_id"]
	return(id_)

def get_bed(line):
	id_ = get_id(line)
	chrom = line[0]
	start = str(int(line[3])-1)
	end = line[4]
	score = "0"
	strand = line[6]
	bed = "\t".join([chrom, start, end, id_, score, strand+"\n"])
	return(bed)

def print_line_to_bed(sp, line):
	bed = get_bed(line)
	feature = line[2]
	output_file = open("./BEDs/{}.{}.bed".format(sp,feature), 'a')
	output_file.write(bed)

def parse_GTF(species):	
	create_dir("./BEDs")
	for sp in species:
		print("Generating bed files for {}".format(sp))
		path=species[sp]["path"]
		os.system("bash ./generate_BEDs.sh {} {}".format(sp, path))
		
		#if "gz" in path:
		#	gtf=gzip.open(path,'rt')
		#else:
		#	gtf=open(path, 'r')
		#for line in gtf:
			#Skip headers and lines other than gene, transcript or exon
		#	if line.startswith("#"): continue
		#	line=line.strip().split("\t")
		#	if not line[2] in ["gene", "transcript", "exon"]: continue
		#	print_line_to_bed(sp, line)

#############################################################################

config=sys.argv[1]
#open("err.log",'w')

#Read config file
species = read_config(config)
names = species.keys()
pairs = itertools.combinations(names, 2)
#print(species, names, pairs)
for p in pairs:
	sp1,sp2=p[0],p[1]
	sp1v,sp2v=species[sp1]["version"],species[sp2]["version"]
	#print(p,sp1, sp1v, sp2, sp2v)
	#map.chain from sp1 to sp2
	chainmaps = retrieve_chainmap_FilesNames(sp1v, sp2v)
	chainmaps_urls = retrieve_chainmap_URLs(sp1v, sp2v, chainmaps)
	#print(chainmaps, chainmaps_urls)
	#Download chainmaps
	download_chainmaps(chainmaps_urls)

	#Parse Input GTF files
	parse_GTF(species)

	#arguments = retrieve_liftover_arguments(species, sp1, sp2, sp1v, sp2v, chainmaps, minMatch)
	#print(arguments)
	#Run liftOver
	#print("liftOver: {} to {}".format(sp1, sp2))
	#subprocess.call(arguments)

