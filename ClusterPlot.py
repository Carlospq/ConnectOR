### IMPORTS
import argparse
import os, sys
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
import numpy as np
import networkx as nx
import pandas as pd


### Options from command line
parser = argparse.ArgumentParser(description = "ConnectOR (v.2.0)"\
                                               "By: Carlos Pulido (carlos.pulidoquetglas@unil.ch)")

parser.add_argument('-gid', '--gene_id',
                    dest = 'gene_id',
                    required = True,
                    help = 'Gene ID of gene to be ploted.')

parser.add_argument('-cid', '--cluster_id',
                    dest = 'cluster_id',
                    required = False,
                    help = 'Cluster ID to be ploted.')

parser.add_argument('-g', '--gene',
                    dest = 'gene_level',
                    action = 'store_true',
                    default = False,
                    help = 'Plot results at gene level (default: exon level)')

parser.add_argument('-p', '--input_path',
                    dest = 'input_path',
                    default = './',
                    help = 'Path to ConnectOR outputs folder (default: ./)')

parser.add_argument('-sp1', '--species1',
                    dest = 'species1',
                    default = '',
                    help = 'Name of species 1.')

parser.add_argument('-sp2', '--species2',
                    dest = 'species2',
                    default = '',
                    help = 'Name of species 2.')

parser.add_argument('-n', '--plot_name',
                    dest = 'plot_name',
                    default = './ClusterPlot.pdf',
                    help = 'File name for pdf to store plot.')

parser.add_argument('-s', '--show',
                    dest = 'show',
                    action = 'store_true',
                    default = False,
                    help = 'Opens a windows with plot in matplotlib.')


### Storing user options
options = parser.parse_args()
gid = options.gene_id
correct_gid = 0
cid = options.cluster_id
gene_level = options.gene_level
input_path = options.input_path
sp1 = options.species1
sp2 = options.species2
plot_name = options.plot_name
show = options.show

feature = "gene" if gene_level else "exon"

if (sp1 and not sp2) or (sp2 and not sp1):
    parser.error("--species1 and --species2 must be specified together.")

if (gid and cid) or (not gid and not cid):
    parser.error("Only one of --gene_id or --cluster_id can be specified (at least one of them must be specified).")


### FUNCTIONS
def geneID_map_to_dict(file_name, dl = "\t"):
    dict_ = {}
    with open(file_name, 'r') as f:
        for line in f:
            line = line.strip().split(dl)
            dict_[line[0]] = {"gene_name": line[1],
                              "gene_type": line[2],
                              "coordinates": line[3],
                              "transcript_IDs": line[4].split("|")}
    return(dict_)


### PLOT
### GeneMaps into dictionary by species
gene_type = {}
sps = []
for input_file in [x for x in os.listdir("./maps")]:
    sp = input_file.split(".")[0]
    sps.append(sp)
    gene_type[sp] = geneID_map_to_dict('maps/%s.geneID_map.txt'%(sp))
    if gid and gid in gene_type[sp]: correct_gid = 1

if not correct_gid:
    raise ValueError( 'Gene ID %s seems not to be in the list of Gene IDs'%(gid) )
    
### Generate clusters
ALL_GENES = []
EDGES = {"exon": [],
         "gene": []}

orthology_path = input_path+'./orthology'
if sp1:
    orthology_files = [x for x in os.listdir(orthology_path) if sp1 in x and sp2 in x]
else:
    orthology_files = [x for x in os.listdir(orthology_path)]

for file in orthology_files:

    if not feature in file: continue

    #Species being analyzed
    org_source = file.split('.')[0].split('to')[0]
    org_destin = file.split('.')[0].split('to')[1]

    #Add all genes to list
    df = pd.read_csv("./maps/%s.geneID_map.txt"%(org_source), sep="\t", header=None,
                     names = ["GeneA_ID", "GeneA_name", "Gene_biotype", "coordinates", "transcripts"])
    ALL_GENES += [(org_source, x, gene_type[org_source][x]["gene_type"]) for x in list(df['GeneA_ID'].unique())]
    ALL_GENES = sorted(set(ALL_GENES))

    #Read pd's with classification
    df = pd.read_csv("./orthology/%s"%(file), sep="\t", header=0,
                     names = ["GeneA_ID", "GeneA_biotype", "GeneB_IDs", "GeneB_biotypes", "GeneB_Classification", "GeneA_number_overlapping_exons", "GeneA_percent_overlapping"])
    df = df[['GeneA_ID','GeneA_biotype', 'GeneB_IDs']]

    #Prepare df's
    df['GeneB_IDs'] = df['GeneB_IDs'].apply(lambda x:[y for y in x.split(',') if y != '.'])

    #Create edges (gene connections/predictions)
    for row in df.itertuples():
        geneA_ID = row.GeneA_ID
        geneA_biotype = row.GeneA_biotype

        for geneB_ID in row.GeneB_IDs:
            geneB_biotype = gene_type[org_destin][geneB_ID]["gene_type"]
            edge = ((org_source, geneA_ID, geneA_biotype), (org_destin, geneB_ID, geneB_biotype))
            EDGES[feature].append(edge)

# ### Create network (clusters)
G_feature = nx.DiGraph()
G_feature.add_nodes_from(ALL_GENES)
G_feature.add_edges_from(EDGES[feature])
G_feature_undirected = G_feature.to_undirected()
components_feature = sorted(list(G_feature_undirected.subgraph(c) for c in nx.connected_components(G_feature_undirected)),key=lambda x:-len(x))
print('%s graph: %d nodes, %d edges, %d connected_components'%(feature, nx.number_of_nodes(G_feature),
                                                               nx.number_of_edges(G_feature), len(components_feature)), file=sys.stderr)

## Plotting component
col_dictionary = {}
palette = plt.get_cmap('Pastel1')
new_palette = ListedColormap(palette(np.arange(len(sps))))
for i in range(len(sps)):
    if sp1:
        if sps[i] not in [sp1, sp2]: continue
    col_dictionary[sps[i]] = new_palette.colors[i]

if gid:
    for c in components_feature:
        for node in c:
            if gid in node:
                component = c
else:
    component = components_feature[cid]

def plot_component(g, file, dictionary, show=show):
    org_color = dictionary
    fig,ax = plt.subplots(figsize=(8,8))
    nodes = g.nodes()
    edges = g.edges()
    edges_set = set(edges)
    bidirectional = [True if (v[1],v[0]) in edges_set else False for v in edges]
    edge_color = ['red' if x else 'gray' for x in bidirectional]
    node_color = [org_color[x[0]] for x in nodes]
    labels = {x:x[1] for x in nodes} 
    pos=nx.spring_layout(g)
    nx.draw_networkx_nodes(g,pos,nodelist=nodes,node_color=node_color)
    nx.draw_networkx_edges(g,pos,edgelist=edges,edge_color=edge_color)
    nx.draw_networkx_labels(g,pos,labels,font_weight='bold')
    plt.legend(handles=[mpatches.Patch(color=v, label=k) for k,v in org_color.items()],
               bbox_to_anchor=(1, 1))
    plt.tight_layout()
    plt.axis('off')
    plt.margins(x=0.3, y=0.3)
    plt.savefig(file, format="pdf")
    if show:
        plt.show()
    
Ge = G_feature.subgraph(component)
plot_component(Ge, plot_name, dictionary=col_dictionary, show=show)

