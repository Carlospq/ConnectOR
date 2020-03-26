import sys

_05dict_ = {}
_1ndict_ = {}
_05map = open(sys.argv[1], 'r')
_1nmap = open(sys.argv[2], 'r')

#read .05 overlap into dictionary
for line in _05map:

	line = line.strip().split("\t")
	t_id = line[0]
	gene_name = line[1]

	_05dict_[t_id] = gene_name
	
#read 1nt overlap into dictionary
for line in _1nmap:

	line = line.strip().split("\t")
	t_id = line[0]
	gene_name = line[1]

	_1ndict_[t_id] = gene_name

for t_id in _05dict_:
	if "novel" in _05dict_[t_id]:
		_05dict_[t_id] = _1ndict_[t_id]

	print "\t".join([t_id, _05dict_[t_id]])