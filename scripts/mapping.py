import sys
import operator

mapp = {}

#
# mapp = { Transcript1: { Exon1: {Gene_name1: 2,
#                                 Gene_name2: 1},
#                         Exon2: {Gene_name1: 1}},
#          Transcript2: { Exon1: {Gene_name3: 1}}}
#

intersections = open(sys.argv[1], 'r')


for line in intersections:
	line = line.strip().split()

	t_id = line[3]
	exon = line[3]+"."+line[4]
	gene_name = line[9]

	if not t_id in mapp:
		mapp[t_id] = {}

	if not exon in mapp[t_id]:
		mapp[t_id][exon] = {}

	if not gene_name in mapp[t_id][exon]:
		mapp[t_id][exon][gene_name] = 1
	else:
		mapp[t_id][exon][gene_name] += 1

counts = {}
for transcript in mapp:

	if not transcript in counts:
		counts[transcript] = {}

	for exon in mapp[transcript]:

		for gene_name in mapp[transcript][exon]:

			if not gene_name in counts[transcript]:
				counts[transcript][gene_name] = mapp[transcript][exon][gene_name]
			else:
				counts[transcript][gene_name] += mapp[transcript][exon][gene_name]

for transcript in counts:

	if '.' in counts[transcript] and len(counts[transcript])>1:
		del counts[transcript]['.'] 

	sorted_x = sorted(counts[transcript].items(), key=operator.itemgetter(1), reverse=True)
	gn = sorted_x[0][0]
	if gn == ".":
		gene_name = ".".join([transcript,"novel"])
	else:
		gene_name = gn

	print "\t".join([transcript, gene_name])
