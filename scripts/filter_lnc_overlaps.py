import sys
import operator

mapp = {}
exons = {}
intersections = open(sys.argv[1], 'r')

#
# mapp = { Transcript1: { Exon1: {Gene_name1: 2,
#                                 Gene_name2: 1},
#                         Exon2: {Gene_name1: 1}},
#          Transcript2: { Exon1: {Gene_name3: 1}}}
#

for line in intersections:

	line = line.strip().split("\t")
	tr_id = line[3]
	g_name = line[9]
	#Exon from sinergia
	exon = line[:6]

	#Initiate dictionary when novel gene found
	if not tr_id in mapp:
		mapp[tr_id] = [line, g_name]
	if not tr_id in exons:
		exons[tr_id] = [exon]

	#Add to dicctionary g_name for each exonic overlap
	if not g_name in mapp[tr_id]:
		mapp[tr_id].append(g_name)
	if not exon in exons[tr_id]:
		exons[tr_id].append(exon)

for tr_id in mapp:
	if len(mapp[tr_id]) == 2:
		if "." in mapp[tr_id]:
			for exon in exons[tr_id]:
				print "\t".join(exon)

