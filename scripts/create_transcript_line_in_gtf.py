import sys

gtf = open(sys.argv[1], 'r')

trs={}
for line in gtf:
	if line.startswith("#"):
		print(line)

	line = line.strip().split("\t")
	arguments = line[-1].split(" ")
	tr_id = arguments[1]
	gen_id = arguments[3]
	chrom = line[0]
	start = int(line[3])
	end = int(line[4])
	strand = line[6]

	if not tr_id in trs:
		trs[tr_id] = [chrom, start, end, strand, gen_id]
	else:
		if start < trs[tr_id][1]:
			trs[tr_id][1] = start
		if end > trs[tr_id][2]:
			trs[tr_id][2] = end

for tr_id in trs:
	chrom = trs[tr_id][0]
	start = str(trs[tr_id][1])
	end = str(trs[tr_id][2])
	strand = trs[tr_id][3]
	gen_id = trs[tr_id][4]

	print(chrom+"\thts\ttranscript\t"+start+"\t"+end+"\t0.000000\t"+strand+"\t.\ttranscript_id "+tr_id+" gene_id "+gen_id)
