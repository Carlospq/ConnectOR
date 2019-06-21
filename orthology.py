import sys

config=open(sys.argv[1], 'r')
species={}
names=[]

#Read config file
for line in config:
	#skip header
	if line.startswith("species_name"): continue

	#generate dictionary with the info for each specie
	line=line.strip().split("\t")
	n,v,p,f=line[0],line[1],line[2],line[3]
	names.append(n)
	species[n]={"path": p, "version": v, "format": f}

#Generate paris to be compared
pairs=[]
for i in range(0, len(names)):
	for j in range(0, len(names)):
		if j <= i: continue
		pairs.append([names[i],names[j]])

for p in pairs:
	sp1,sp2=p[0],p[1]

	#map.chain from sp1 to sp2
	map_chain=sp1+"To"+sp2+".over.chain.gz"
	url_map_chain='http://hgdownload.cse.ucsc.edu/goldenPath/{}/liftOver/{}'.format(sp1, map_chain)


#olfFile = BED/GTF/GFF file with annotation to liftOver
#oldFile=species[$sp1]

#newFile=$4


#unMapped=$5


#minMatch=0.95

#liftOver -minMatch $minMatch $oldFile $map_chain $newFile $unMapped

