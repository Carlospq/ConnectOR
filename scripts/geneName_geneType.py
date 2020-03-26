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

dic_ = {}
mapp = {}
tgmap = sys.argv[1]
gtf = sys.argv[2]

#GENCODE (or another reference gtf) with "gene_name" flag
gtf = open(gtf, 'r')
for line in gtf:

	if line.startswith("#"): continue
	line = line.strip().split("\t")
	if not "gene" == line[2]: continue
	arguments = arguments_dic(line[-1].split(" "))

	gene_name = arguments["gene_name"]
	gene_type = arguments["gene_type"]

	if not gene_name in dic_:
		dic_[gene_name] = [gene_type]
	else:
		dic_[gene_name].append(gene_type)

#Read mapp (transcript_id - gene_name)
for line in open(tgmap, 'r'):
	line = line.strip()
	line = re.split(' |\t',line)
	if line[1] in dic_:
		mapp[line[1]] = dic_[line[1]]
	else:
		mapp[line[1]] = ["NOVEL"]

for gene_name in mapp:
	types = ",".join(set(sorted(mapp[gene_name])))
	print "\t".join([gene_name,types])