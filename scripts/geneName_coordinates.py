import sys, re

##### function to parse GTF arguments into a dictionary #####
def arguments_dic(line):
	dic = {}
	n = 0
	for e in line:
		if n == 0:
			key = e
			n = 1
			continue
		if n == 1:
			dic[key] = e.replace("\"","").replace(";","")
			n = 0
	return dic
#############################################################

mapp = {}
dic_ = {}

tgmap = sys.argv[1]
gtf = sys.argv[2]

#Read mapp (transcript_id - gene_name)
for line in open(tgmap, 'r'):
	line = line.strip()
	line = re.split(' |\t',line)
	mapp[line[0]] = line[1]

#Read gtf file to create gene_name - coordinates mapp
gtf = open(gtf, 'r')
for line in gtf:

	if line.startswith("#"): continue
	line = line.strip().split("\t")
	if not "transcript" == line[2]: continue
	arguments = arguments_dic(line[-1].split(" "))

	chrom = line[0]
	start = int(line[3])
	end = int(line[4])
	strand = line[6]

	transcript_id = arguments["transcript_id"]
	gene_name = mapp[transcript_id]

	if not gene_name in dic_:
		dic_[gene_name] = [chrom, start, end, strand]
	else:
		if start < dic_[gene_name][1]:
			dic_[gene_name][1] = start
		if end > dic_[gene_name][2]:
			dic_[gene_name][2] = end

#Print Gene coordinates in BED format
for gid in dic_:
	chrom = dic_[gid][0]
	start = str(dic_[gid][1])
	end = str(dic_[gid][2])
	strand = dic_[gid][3]
	print "\t".join([gid, chrom, start, end, strand])

