### IMPORTS
import re
import subprocess
import json
import wget
import gzip
import itertools
import random
import pickle
import os, sys
from os import path,listdir
from io import StringIO
from tqdm import tqdm
from collections import Counter
import fileinput
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

#Only for jupyter
#from IPython.display import IFrame

### Dictionaris
with open('dictionaries.json') as f:
    dictionaries = json.load(f)
    
def check_folder(folder_name):
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
		

### VARIABLES
# minMatch liftOver required to mapp to new region
try:
    minMatch=sys.argv[1]
    if int(minMatch) < 0 or int(minMatch) > 99:
        sys.exit('minMatch for LiftOver must be between 0 and 99')
except IndexError:
    minMatch=50

### FUNCTIONS
def assign_geneType(Type):
    type_dict = {"none": ["none"],
                 "pc": ["protein_coding", "pc"],
                 "ncRNA": ["NOVEL", "stringtie", "3prime_overlapping_ncRNA", "antisense", "bidirectional_promoter_lncRNA", "lincRNA", "ncRNA",
                           "non_coding", "processed_transcript", "sense_intronic", "sense_overlapping", "lncRNA", "snoRNA", "snRNA", "sRNA", 
                           "misc_RNA", "rRNA", "scaRNA", "miRNA", "TEC", "scRNA", "macro_lncRNA",
                           "pseudogene", "transcribed_processed_pseudogene", "translated_processed_pseudogene", "polymorphic_pseudogene", "processed_pseudogene", "rRNA_pseudogene",
                           "transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene", "translated_unprocessed_pseudogene"],
                 "other": ["other", "IG_V_pseudogene", "IG_C_pseudogene", "IG_V_gene", "IG_C_gene", "IG_J_gene", "TR_J_pseudogene", "IG_J_pseudogene", "IG_pseudogene", "IG_D_pseudogene", 
                           "ribozyme", "vault_RNA", "vaultRNA", "TR_C_gene", "TR_J_gene", "TR_V_gene", "TR_V_pseudogene", "TR_D_gene", "IG_D_gene", "Mt_tRNA", "Mt_rRNA", "IG_LV_gene"]}
    biotype = [t for t in type_dict if Type in type_dict[t]]
    return(biotype[0])

def read_config(file_name):
    config_df =  pd.read_csv(file_name, sep='\t', na_filter= False)
    return(config_df)


def check_folder(folder_name):
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

        
def check_file(file_name_list):
    for file_name in file_name_list:
        if not os.path.isfile(file_name):
            print("'%s' doesn not exist. Please check config_file. Exiting..."%(file_name))
            sys.exit(1)

            
def download_default_files(file_url, folder_name):
    check_folder(folder_name)
    file_name = file_url.split("/")[-1]
    if not os.path.isfile("/".join([folder_name, file_name])):
        wget.download(file_url, folder_name)
        print('\t%s downloaded succesfully..'%(file_name))
    else:
        print('\t%s already exists.. Skipping'%(file_name))

        
def arguments_dic(line):
    args = line[-1].split(";")
    dic = {}
    for e in args:
        if e == ";": continue
        if e == "": continue
        if e[0] == " ": e = e[1:]
        key = e.split(" ")[0].replace("\"","")
        dic[key] = e.split(" ")[1].replace("\"","")
    return dic


def generate_maps(gtf, sp):
    
    if gtf.endswith(".gz"):
        f = gzip.open(gtf, 'rb')
        compressed = True
    else:
        f = open(gtf, 'r')
        compressed = False
        
    transcripts = {}
    genes = {}

    for line in f:
        if compressed: line = str(line, 'utf-8')
        if line.startswith("#"): continue
        line = line.strip().split("\t")
        arguments = arguments_dic(line)

        if "gene_name" in arguments:
            gene_name = arguments["gene_name"]
        else:
            gene_name = arguments["gene_id"]
            
        gene_biotype = biotype(["gene_type", "gene_biotype"], arguments)
        gene_biotype = assign_geneType(gene_biotype)
        if line[2] == "transcript":
            transcripts[arguments["transcript_id"]] = {"gene_id": arguments["gene_id"],
                                                       "gene_type": gene_biotype,
                                                       "gene_name": gene_name}
            if arguments["gene_id"] in genes:
                if genes[arguments["gene_id"]]["gene_type"] == "protein_coding":
                    gene_biotype = "protein_coding"
            genes[arguments["gene_id"]] = {"gene_type": gene_biotype,
                                           "gene_name": gene_name}
            
    fo1 = open("maps/"+sp+".transcriptID_geneID_map.txt", 'w')
    for t in transcripts:
        line = "\t".join([t, transcripts[t]["gene_id"], transcripts[t]["gene_name"], transcripts[t]["gene_type"]])
        fo1.write(line+"\n")
    fo1.close()
    
    if len(genes) > 0:
        fo2 = open("maps/"+sp+".geneID_geneName_geneType_map.txt", 'w')
        for g in genes:
            line = "\t".join([g, genes[g]["gene_name"], genes[g]["gene_type"]])
            fo2.write(line+"\n")
        fo2.close()

        
def biotype(keys, arguments):
    biotype=''
    while not biotype:
        for k in keys:
            try:
                biotype = arguments[k]
                return(biotype)
            except KeyError:
                pass
        if not biotype: biotype = "NOVEL"
        return(biotype)

    
def generate_beds(gtf, sp):
    
    if gtf.endswith(".gz"):
        f = gzip.open(gtf, 'rb')
        compressed = True
    else:
        f = open(gtf_file, 'r')
        compressed = False
        
    transcripts = {}
    genes = {}

    output_exons = open("./BEDs/{}.exons.bed".format(sp),"w")
    for line in f:
        if compressed: line = str(line, 'utf-8')
        if line.startswith("#"): continue
        line = line.strip().split("\t")

        if line[2] != "exon": continue
        arguments = arguments_dic(line)

        g_id = arguments["gene_id"]
        chrom = line[0] if line[0].startswith("chr") else "chr"+line[0]
        start = str(int(line[3])-1)
        end = line[4]
        strand = line[6]
        if int(start) < 0: start = "0"
        if not g_id in genes:
            genes[g_id] = {"chrom": chrom, 
                           "start": int(start),
                           "end": int(end),
                           "strand": strand,
                           "gene_name": ""}
        else:
            if int(start) < genes[g_id]["start"]:
                genes[g_id]["start"] = int(start)
            if int(end) > genes[g_id]["end"]:
                genes[g_id]["end"] = int(end)

        if "gene_name" in arguments:
            genes[g_id]["gene_name"] = arguments["gene_name"]
        else:
            genes[g_id]["gene_name"] = arguments["gene_id"]

        exon_bed_line = "\t".join([chrom, start, end, genes[g_id]["gene_name"], '0', strand])+"\n"
        output_exons.write(exon_bed_line)
    output_exons.close()
    
    output_genes = open("./BEDs/{}.genes.bed".format(sp),"w")
    for gene in genes:
        d = genes[gene]
        chrom = d["chrom"]
        start = str(d["start"])
        end = str(d["end"])
        strand = d["strand"]
        gene_name = d["gene_name"]

        gene_bed_line = "\t".join([chrom, start, end, gene_name, '0', strand])+"\n"
        output_genes.write(gene_bed_line)
    output_genes.close()
    
    
def bed_sort(sp_v):
    files = ["./BEDs/{}.exons.bed".format(sp_v), "./BEDs/{}.genes.bed".format(sp_v)]
    for file_name in files:
        print("\r\t"+file_name+"... sorting", end="")
        call = "sort -u -k1,1 -k2,2n -o '%s' '%s'"%(file_name,file_name)    
        subprocess.call(call, shell=True)
        print("\r\t"+file_name+"... sorted ")
        
        
def gene_map_to_dict(file_name):
    d = {}
    finput = fileinput.FileInput(files=file_name)
    for line in finput:
        line = line.strip().split("\t")
        d[line[1]] = {"gene_id": line[0],
                      "gene_name": line[1],
                      "gene_type": line[2]}
    finput.close()
    return(d)


def transcript_map_to_dict(file_name, dl = "\t"):
    dict_ = {}
    with open(file_name, 'r') as f:
        for line in f:
            line = line.strip().split(dl)
            dict_[line[2]] = {"transcript_ID": line[0],
                              "gene_id": line[1],
                              "gene_type": line[3]}
    return(dict_)


def parse_orthologs(line):
    geneM = line[0]
    atype = line[1]
    ortho = line[2].split(",")
    nexon = line[3].split(",")
    pcent = line[4].split(",")
    btype = line[5].split(",")

    return(geneM, atype, ortho, nexon, pcent, btype)


def count_classes(btype):
    tmp   = []
    none  = ["none"]
    pc    = ["protein_coding", "pc"]
    ncRNA = ["NOVEL", "stringtie", "3prime_overlapping_ncRNA", "antisense", "bidirectional_promoter_lncRNA", "lincRNA",
             "non_coding", "processed_transcript", "sense_intronic", "sense_overlapping", "lncRNA", "snoRNA", "snRNA", "sRNA", 
             "misc_RNA", "rRNA", "scaRNA", "miRNA", "TEC", "scRNA", "macro_lncRNA",
             "pseudogene", "transcribed_processed_pseudogene", "translated_processed_pseudogene", "polymorphic_pseudogene", "processed_pseudogene", "rRNA_pseudogene",
             "transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene", "translated_unprocessed_pseudogene"]
    other = ["IG_V_pseudogene", "IG_C_pseudogene", "IG_V_gene", "IG_C_gene", "IG_J_gene", "TR_J_pseudogene", "IG_J_pseudogene", "IG_pseudogene", "IG_D_pseudogene", 
             "ribozyme", "vault_RNA", "TR_C_gene", "TR_J_gene", "TR_V_gene", "TR_V_pseudogene", "TR_D_gene", "IG_D_gene", "Mt_tRNA", "Mt_rRNA", "IG_LV_gene"]
    for type_ in btype:
        if type_ in none:
            tmp.append("none")
        elif type_ in ncRNA:
            tmp.append("lncRNA")
        elif type_ in pc:
            tmp.append("pc")
        else:
            tmp.append("other")
    class_ = [tmp.count("none"), tmp.count("lncRNA"), tmp.count("pc"), tmp.count("other")]

    return(class_)


def classification(classes):
    c = "other"
    #unique cases
    #  none              lncRNA            pc                other
    if classes[0]==1 and classes[1]==0 and classes[2]==0 and classes[3]==0:
        c = "none"
    if classes[0]>=0 and classes[1]==1 and classes[2]==0 and classes[3]==0:
        c = "lncRNA"
    if classes[0]>=0 and classes[1]==0 and classes[2]==1 and classes[3]==0:
        c = "pc"
    if classes[0]>=0 and classes[1]==0 and classes[2]==0 and classes[3]==1:
        c = "other"

    #multiple cases
    #  none              lncRNA            pc                other
    if classes[0]>1  and classes[1]==0 and classes[2]==0 and classes[3]==0:
        c = "nones"
    if classes[0]>=0 and classes[1]>1  and classes[2]==0 and classes[3]==0:
        c = "lncRNAs"
    if classes[0]>=0 and classes[1]==0 and classes[2]>1  and classes[3]==0:
        c = "pcs"
    if classes[0]>=0 and classes[1]==0 and classes[2]==0 and classes[3]>1:
        c = "others"

    #dual cases
    #  none              lncRNA            pc                other
    if classes[0]>=0 and classes[1]>=1 and classes[2]>=1 and classes[3]>=0:
        c = "lncRNA_PC"
    if classes[0]>=0 and classes[1]>=1 and classes[2]==0 and classes[3]>=1:
        c = "lncRNA_other"
    if classes[0]>=0 and classes[1]==0 and classes[2]>=1 and classes[3]>=1:
        c = "PC_other"

    return(c)



def change_default(feature, value, config_df, i):
    default_values = config_df.loc[i]["default"].split("|")
    if feature == "gtf": 
        n = 0
    else:
        n = 1
    default_values[n] = value
    return("|".join(default_values))


def assign_class(df, col_check1 = "enames", col_check2 = "gnames", col_assign = "class"):
    for index, row in df.iterrows():

        set1 = set(row[col_check1].split(";"))
        set2 = set(row[col_check2].split(";"))

        for _set_ in [set1, set2]:
            if "." in _set_:
                _set_ = _set_.remove(".") if len(_set_) > 1 else _set_

        #Defining Class
        #Sets are equal
        if (set1 == set2):
            #Equal one to one
            if (len(set1)==1) & (len(set2)==1):
                #Sets are empty (".")
                if ("." in set1) & ("." in set2):
                    row[col_assign] = "class6"
                #Sets have same gene_name
                else:
                    row[col_assign] = "class1"
            #Equal many to many
            elif (len(set1)>1) & (len(set2)>1):
                row[col_assign] = "class2"
        #Sets are different
        elif (set1 != set2):
            #One set is == "."
            if ((set1 == {"."}) & (set2 != {"."})) | ((set1 != {"."}) & (set2 == {"."})):
                row[col_assign] = "class5"
            #eclass is subset of gclass
            if (set1 <= set2) & (not set1 >= set2):
                row[col_assign] = "class3"
            #eclass is superset of gclas (*should not happen)
            elif (set1 >= set2) & (not set1 <= set2):
                row[col_assign] = "class4"

    return(df)


# SUBCLASS
subclass = {"lncRNA lncRNA":             "subclass1",
            "lncRNA lncRNAs":            "subclass2",
            "lncRNAs lncRNAs":           "subclass3",
            "lncRNA_other lncRNA_other": "subclass4",
            
            "pc pc":                     "subclass5",
            "pc pcs":                    "subclass6",
            "pcs pcs":                   "subclass7",
            "pc pc_other":               "subclass8",
            "pc_other pc_other":         "subclass9",

            "lncRNA lncRNA_pc":          "subclass10",
            "lncRNAs lncRNA_pc":         "subclass11",
            "lncRNA_pc lncRNA_pc":       "subclass12",
            
            "pc lncRNA_pc":              "subclass13",
            "pcs lncRNA_pc":             "subclass14",
            "pc_other lncRNA_pc":        "subclass15",
            
            "other other":               "subclass16",
            "others others":             "subclass17",
            "other lncRNA_other":        "subclass18",
            "other pc_other":            "subclass19",
            
            "none none":                 "subclass20",
            "none lncRNA":               "subclass21",
            "none lncRNA_pc":            "subclass22",
            "none lncRNAs":              "subclass23",
            "none pc":                   "subclass24",
            "none pcs":                  "subclass25"}


def assign_subclass(df, sc_dict=subclass, col_check1 = "eclass", col_check2 = "gclass", col_assign = "subclass"):
    df[col_assign] = df[col_check1] + " " + df[col_check2]
    #df = df.replace({col_assign: sc_dict})
    return(df)
	
def get_cluster(sp, gene_name, g, bidirectionality = False, sp2 = False):
    #Returns Cluster of specific gene and all its connections
    #bidirectionality = True: returns only nodes instead of edges and counts of ocurrences for each node
    #If sp2 != False: 
    biotype = gene_type[sp][gene_name]["gene_type"]
    node = (sp, gene_name, biotype)
    
    s = g.subgraph(nx.shortest_path(g.to_undirected(), node))
    edges = s.edges()    
    edges_set = set(edges)

    nodes = []
    links = []
    
    for edge in edges:
        sps_in_edge = [node[0] for node in edge]
        if sp2:
            if not sp in sps_in_edge or not sp2 in sps_in_edge: continue
        is_bidirectional = 1 if (edge[1],edge[0]) in edges_set else 0
        l1 = list(edge)
        #1=edge is bidirectional; 2=edge is unidirectional
        l1.append(is_bidirectional if is_bidirectional else 2)
        links.append(tuple(l1))
        
        n1 = list(edge[0])
        n2 = list(edge[1])
        #1=gene has bidirectional prediction; 2=gene has unidirectional prediction; 3=gene has no prediction
        if is_bidirectional:
            n1.append(1)
            nodes.append(tuple(n1))
        else:
            n1.append(2)
            n2.append(3)           
            nodes.append(tuple(n1))
            nodes.append(tuple(n2))
    
    if bidirectionality:
        #bidirectionality = True: returns only nodes instead of edges and counts of ocurrences for each node
        return(Counter(nodes))
    else:
        return(links)
		
def generate_components_stats(components, G):
    x_to_comp_idx = {}
    comp_idx_x = {}
    for idx,c in enumerate(components):

        s = G.subgraph(c)
        edges = s.edges()    
        edges_set = set(edges)
        nodes = []
        links = []
        for gene in c:
            #gene[1]=GeneName; idx=cluster_id; gene=tuple(sp, GeneName, GeneBiotype) 
            x_to_comp_idx[gene[1]] = {"idx": idx, "node": gene}

        for edge in edges:

            sps_in_edge = [node[0] for node in edge]
            is_bidirectional = 1 if (edge[1],edge[0]) in edges_set else 0
            l1 = list(edge)
            #1=edge is bidirectional; 2=edge is unidirectional
            l1.append(is_bidirectional if is_bidirectional else 2)
            links.append(tuple(l1))

            n1 = list(edge[0])
            n2 = list(edge[1])
            #1=gene has bidirectional prediction; 2=gene has unidirectional prediction; 3=gene has no prediction
            if is_bidirectional:
                n1.append(1)
                nodes.append(tuple(n1))
            else:
                n1.append(2)
                n2.append(3)           
                nodes.append(tuple(n1))
                nodes.append(tuple(n2))

        #bidirectionality = True: returns only nodes instead of edges and counts of ocurrences for each node
        comp_idx_x[idx] = {"bidirectional": Counter(nodes),
                           "links": links}
    return(x_to_comp_idx, comp_idx_x)
	

def plot_component(g, d):
    org_color = d
    fig,ax = plt.subplots(figsize=(8,8))
    nodes = g.nodes()
    edges = g.edges()    
    edges_set = set(edges)
    bidirectional = [True if (v[1],v[0]) in edges_set else False for v in edges]
    edge_color = ['red' if x else 'gray' for x in bidirectional]
    node_color = [org_color[x[0]] for x in nodes]
    labels = {x:x[1] for x in nodes} ############### cambiar x[1] por tu_dictionary[x[1]]!
    pos=nx.spring_layout(g)    
    nx.draw_networkx_nodes(g,pos,nodelist=nodes,node_color=node_color)
    nx.draw_networkx_edges(g,pos,edgelist=edges,edge_color=edge_color)
    nx.draw_networkx_labels(g,pos,labels,font_weight='bold')
    plt.legend(handles=[mpatches.Patch(color=v, label=k) for k,v in org_color.items()],
               bbox_to_anchor=(1, 1))
    plt.tight_layout()
    plt.axis('off')
    plt.show()
	
def check_connections(links, sp1, sp2):

    edges = [tuple([edge[0], edge[1]]) for edge in links]
    nodes = []
    for edge in edges:
            
        sps_in_edge = [edge[0][0], edge[1][0]]
        if not sp1 in sps_in_edge or not sp2 in sps_in_edge: continue
        is_bidirectional = 1 if (edge[1],edge[0]) in edges else 0

        n1 = list(edge[0])
        n2 = list(edge[1])
        #1=gene has bidirectional prediction; 2=gene has unidirectional prediction; 3=gene has no prediction
        if is_bidirectional:
            n1.append(1)
            nodes.append(tuple(n1))
        else:
            n1.append(2)
            n2.append(3)           
            nodes.append(tuple(n1))
            nodes.append(tuple(n2))

    return(Counter(nodes))
    
    
def get_clusters_stats(components, G, level, output_path="./counts/", species_names=False):
    if species_names:
        species = species_names
    else:
        species = sps.keys()
        
    #Biotypes accepted for ncRNA clusters
    ncrna_bt = ["NOVEL", "stringtie", "3prime_overlapping_ncRNA", "antisense", "bidirectional_promoter_lncRNA", "lincRNA", "ncRNA",
                "non_coding", "processed_transcript", "sense_intronic", "sense_overlapping", "lncRNA", "snoRNA", "snRNA", "sRNA", 
                "misc_RNA", "rRNA", "scaRNA", "miRNA", "TEC", "scRNA", "macro_lncRNA",
                "pseudogene", "transcribed_processed_pseudogene", "translated_processed_pseudogene", "polymorphic_pseudogene", "processed_pseudogene", "rRNA_pseudogene",
                "transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene", "translated_unprocessed_pseudogene"]
    
    #Check input/output names
    if not path.exists(output_path):
        os.makedirs(output_path)
    
    if species_names:
        o_name = "".join([species_names[0],"_", species_names[1], "_"])
    else: 
        o_name = ""
        
    output_cluster_stats = path.join(output_path, '%scluster_stats_%s.csv'%(o_name, level))
    output_cluster_genes = path.join(output_path, '%scluster_genes_%s.csv'%(o_name, level))    
    
    org_all_count = list(species)
    
    #Total genes
    total_genes = len(G.nodes())
    
    #Headers for output files
    header_stats = ['Cluster ID','Nodes','Number of species', 'Bidirectionality', 'Cluster type', 'Biotypes']
    header_genes = ['Gene Name','Species', 'Biotype', 'Cluster ID', 'Cluster Biotype', 'Cluster type', 'Orthologues', 'in_degree','out_degree']

    for sp in species:
        header_stats += ['has_'+sp, 'count_'+sp, 'Gene Names '+sp]
        header_genes.insert(6, 'Gene to '+sp)
    
    o1 = open(output_cluster_stats,'w')
    o1.write(','.join(header_stats)+'\n')
    o2 = open(output_cluster_genes,'w')
    o2.write(','.join(header_genes)+'\n')
    
    #Iterate over clusters
    for idx,c in enumerate(components):
        g = G.subgraph(c)
        nodes = g.nodes()
        edges = g.edges()    
        edges_set = set(edges)
        print("\r%.2f%s Cluster %s of %s (%s genes in cluster)"%(int(idx)/len(components)*100, "%", str(idx), str(len(components)), str(len(g.nodes()))) , end= "\r", file=sys.stderr)
        #Split in nodes and edges


        #Bidirectionality in the cluster
        is_bidirectional = [1 if (v[1],v[0]) in edges_set else 0 for v in edges]
        bidirectional = 0 if len(edges) == 0 else sum(is_bidirectional)/len(edges)

        #Biotypes in cluster
        biotypes = [n[2] for n in nodes]
        for bt in biotypes:
            if bt not in ncrna_bt:
                biotypes = "other"
                break
            else:
                biotypes = "ncRNA"
        
        #Nodes per species in cluster
        org_comp_count_ugly = Counter([x[0] for x in nodes])
        org_comp_count = {k:v for k,v in org_comp_count_ugly.items()}
        for sp in species:
            if sp not in org_comp_count:
                org_comp_count[sp] = 0
        n_species = [1 for sp in org_comp_count if org_comp_count[sp] > 0]
        n_species = sum(n_species)
        
        #CLuster data
        data = {'Cluster ID': str(idx),
                'Nodes': str(len(nodes)),
                'Number of species': str(n_species),
                'Bidirectionality': '%.2f'%(bidirectional),
                'Biotypes': biotypes}
        
        for sp in species:
            data['has_'+sp] = str(int(org_comp_count[sp]>0))
            data['count_'+sp] = str(org_comp_count[sp])
            data['per_cluster_'+sp] = '%.2f'%(100*org_comp_count[sp]/len(nodes))
            data['per_total_'+sp] = '%.2f'%(100*org_comp_count[sp]/total_genes)
            data['Gene Names '+sp] = "|".join([n[1] for n in nodes if n[0]==sp])
                                              
        #Cluster type
        ctype=''
        has_list = [int(data[x]) for x in data if x.startswith("has")]
        
        count_list = [int(data[x]) for x in data if x.startswith("count")]
        count_list_c1 = 1 if 1 in Counter(count_list) else 0
        count_list_c2 = 0
        
        for c in Counter(count_list):
            if c == 1: continue
            if Counter(count_list)[c] >= 1: count_list_c2 = 1
            
        if len(nodes) <= sum(has_list) and bidirectional == 1: ctype="One to one"
        if len(nodes) <= sum(has_list) and bidirectional <  1: ctype="One to half"
        if len(nodes) == 1: ctype="One to none"
                
        if len(nodes) > sum(has_list) and count_list_c1 and count_list_c2: ctype="One to many"
        if len(nodes) > sum(has_list) and not count_list_c1 and count_list_c2: ctype="Many to many"
        data['Cluster type'] = ctype
        
        o1.write(','.join(data[x] for x in header_stats)+'\n')
       
        #Write info for each node (gene) to output file 2
        comp_idx_x = comp_idx_exon if level == "exon" else comp_idx_gene
        x_to_comp_idx = exon_to_comp_idx if level == "exon" else gene_to_comp_idx
        i=0
        for node in nodes:
            i+=1
            print("\r%.2f%s Cluster %s of %s (%s/%s genes in cluster)             "%(int(idx)/len(components)*100, "%", str(idx), str(len(components)), i, str(len(g.nodes()))) , end= "\r", file=sys.stderr)
            predictions = "|".join([";".join([e[1][1], e[1][0],e[1][2]]) for e in edges if node == e[0]])
            biotype = node[2]
            data2 = {'Gene Name': node[1],
                     'Species': node[0],
                     'Biotype':biotype,
                     'Cluster ID': str(idx),
                     'Cluster Biotype': biotypes,
                     'Cluster type': data['Cluster type'],
                     'Orthologues': predictions,
                     'in_degree': str(g.in_degree(node)),
                     'out_degree': str(g.out_degree(node))}
        
            #Get connections in cluster only for sp1 and sp2
            for sp in species:
                sc = ""
                #sp == sp1: don't analyze
                if sp == node[0]:
                    ort_type = ''
                #sp != sp1: analyze
                else:
                    ort_type = ''
                    in_d  = int(data2['in_degree'])
                    out_d = int(data2['out_degree'])

                    comp_idx_x = comp_idx_exon if level == "exon" else comp_idx_gene
                    x_to_comp_idx = exon_to_comp_idx if level == "exon" else gene_to_comp_idx
                    c = comp_idx_x[x_to_comp_idx[node[1]]["idx"]]["links"]
                    c = check_connections(c, node[0], sp)
                    sc = [s for s in c if s[1] == node[1]]
                    if not sc:
                        n_type = -1
                    else:
                        n_type = sc[0][3] if len(sc)==1 else min([s[3] for s in sc])

                    #Gene with same (one) in and out prediction
                    if in_d == 1 and out_d == 1 and n_type == 1:
                        ort_type = "One to one"
                    #Gene with no in_prediction but out_predictions
                    elif in_d == 0 and out_d == 1:
                        ort_type = "One to half"
                    #Gene with in_prediction but no out_predictions
                    elif in_d == 1 and out_d == 0:
                        ort_type = "One to half"
                    #Gene with (> one) different in_prediction and out_predictions
                    elif in_d >= 0 and out_d >= 0 and n_type >= 0:
                        ort_type = "One to many"            
                    #Gene with no in or out prediction
                    elif in_d == 0 and out_d == 0:
                        ort_type = "One to none" 

                data2["Gene to "+sp] = ort_type

            o2.write(','.join(data2[x] for x in header_genes)+'\n')
    
    o1.close()
    print('Wrote',output_cluster_stats,file=sys.stderr)
    o2.close()
    print('Wrote',output_cluster_genes,file=sys.stderr)

		
### Main
## Dictionaris
with open('dictionaries.json') as f:
  dictionaries = json.load(f)

## Read config file
config_df = read_config("./config")
config_df.set_index("species", inplace = True)
config_df["default"] = "False|False"
for i in config_df.index:
    if not config_df.loc[i]["annotation"]:
        default_gtf = dictionaries["gtfs_ensembl_r98"][config_df.loc[i]["assembly_version"].lower()].split("/")[-1]
        config_df.at[i, 'default'] = change_default("gtf", "True", config_df, i)
        config_df.at[i, 'annotation'] = default_gtf
    if not config_df.loc[i]["chainmap"]:
        chainmaps = []
        config_df.at[i, 'default'] = change_default("chainmap", "True", config_df, i)
        for j in config_df.index:
            if i!=j:        
                default_chainmap_path = dictionaries["chain_maps"][config_df.loc[i]["assembly_version"].lower()][config_df.loc[j]["assembly_version"].lower()]
                default_chainmap_name = "chainmaps/"+default_chainmap_path.split("/")[-1]
                chainmaps.append(default_chainmap_name)
        config_df.at[i, 'chainmap'] = ",".join(chainmaps)

print("Checking Config_file...")
for i in config_df.index:
    if not config_df.loc[i]["assembly_version"].lower() in dictionaries["chain_maps"]:
        print("'%s' is not a valid assembly_version. Please use one of the following values for default analysis: %s"%(config_df.loc[i]["assembly_version"], ",".join(dictionaries["chain_maps"])))
        sys.exit(1)
    
    defaults = config_df.loc[i]["default"].split("|")
    #Check GTFs and chainmaps
    for feature,default in zip(["annotation","chainmap"],defaults):
        sp_vi = config_df.loc[i]["assembly_version"].lower()
        if default == "True":
            if feature == "GTFs": download_default_files(dictionaries["gtfs_ensembl_r98"][sp_vi], "GTFs")
            if feature == "chainmap":
                for j in config_df.index:
                    if i == j: continue
                    sp_vj = config_df.loc[j]["assembly_version"].lower()
                    url=dictionaries["chain_maps"][sp_vi][sp_vj]
                    download_default_files(url, "chainmaps")
        else:
            if feature == "annotation": check_file([config_df.loc[i][feature]])
            if feature == "chainmap":   check_file(config_df.loc[i][feature].split(","))
print("Config_file is correct..")

## Download GTF files (if no gtf provided in config)
print('Downloading GTFs')
for i in config_df.index:
    if config_df.loc[i]["default"].split("|")[0] == "False": 
        print("\tNo default GTF needed for %s.. Skipping"%(i))
        continue    
    sp_v = config_df.loc[i]["assembly_version"].lower()
    download_default_files(dictionaries["gtfs_ensembl_r98"][sp_v], "GTFs")        
	
## Download chainmaps
print('Downloading default chainmaps...')
for i in config_df.index:
    if config_df.loc[i]["default"].split("|")[1] == "True":
        sp_vi = config_df.loc[i]["assembly_version"].lower()
        for j in config_df.index:
            if i == j: continue
            sp_vj = config_df.loc[j]["assembly_version"].lower()
            url=dictionaries["chain_maps"][sp_vi][sp_vj]
            download_default_files(url, "chainmaps")
    else:
        print("\tNo default chainmaps needed for %s.. Skipping"%(i))
		
## Generate maps (transcriptID-geneID & geneID-geneName-geneType)
check_folder("maps")
print("Generating transcriptID-geneID & geneID-geneName-geneType maps...")
for i in config_df.index:
    gtf_file = config_df.loc[i]["annotation"]
    sp_v = config_df.loc[i]["assembly_version"].lower()
    if config_df.loc[i]["default"].split("|")[0] == "True":
        gtf_file = "GTFs/"+gtf_file
    print("\r\tMaps for %s... generating"%(sp_v), end="")
    generate_maps(gtf_file, sp_v)
    print("\r\tMaps for %s... generated "%(sp_v))

	
## Generate BED with genes/exons from GTF files
check_folder("BEDs")
print("Generating BED files for exons and genes...")
for i in config_df.index:
    gtf_file = config_df.loc[i]["annotation"]
    sp_v = config_df.loc[i]["assembly_version"].lower()
    if config_df.loc[i]["default"].split("|")[0] == "True":
        gtf_file = "GTFs/"+gtf_file
    print("\tBEDs for %s... generating"%(sp_v), end="")
    generate_beds(gtf_file, sp_v)
    print("\r\tBEDs for %s... generated "%(sp_v))
	
print("Sorting BED files...")
for i in config_df.index:
    sp_v = config_df.loc[i]["assembly_version"].lower()
    bed_sort(sp_v)
	
## LifOver exons/genes
features = ["exons", "genes"]
check_folder("liftovers")
print("LiftOver...")
for i in config_df.index:
    sp_vi = config_df.loc[i]["assembly_version"].lower()
    chainmaps = config_df.loc[i]["chainmap"].split(",")
    n=0
    for j in config_df.index:
        if i == j: continue
        sp_vj = config_df.loc[j]["assembly_version"].lower()
        map_chain = chainmaps[n] if config_df.loc[i]["default"] else chainmaps[n]
        n+=1
        #liftOver oldFile map.chain newFile unMapped
        for feature in features:
            print("\r\t{} {} to {}... mapping".format(i, feature, j))
            oldFile = "BEDs/{}.{}.bed".format(sp_vi, feature)
            newFile = "liftovers/{}to{}.{}.liftover".format(sp_vi, sp_vj, feature)
            unMapped= "liftovers/{}to{}.{}.unmapped".format(sp_vi, sp_vj, feature)
            os.system("./scripts/liftOver -minMatch=0.{} {} {} {} {}".format(minMatch, oldFile, map_chain, newFile, unMapped))
            print("\r\t{} {} to {}... done         ".format(i, feature, j))

## Intersect LiftOvers
check_folder("overlaps")
print("Intersecting LiftOver...")
for i in config_df.index:
    for j in config_df.index:
        if i == j: continue
        for f in ["exons", "genes"]:
            print("\r\t{} {} to {}... intersecting".format(i, f, j), end="")
            sp1=config_df.loc[i]["assembly_version"].lower()
            sp2=config_df.loc[j]["assembly_version"].lower()
            lifover_input = 'liftovers/%sto%s.%s.liftover'%(sp1, sp2, f)
            bed_input = 'BEDs/%s.%s.bed'%(sp2, f)
            output = 'overlaps/%sto%s.%s.overlap'%(sp1, sp2, f)
            call = 'intersectBed -wao -s -a %s -b %s > %s'%(lifover_input, bed_input, output)
            subprocess.call(call, shell=True, executable='/bin/bash')
            print("\r\t{} {} to {}... done         ".format(i, f, j))

## Parse Overlap
check_folder("orthology")
print("Predicting orthologues...")

for i in config_df.index:
    for j in config_df.index:
        
        if i == j: continue
        
        sp1=config_df.loc[i]["assembly_version"].lower()
        sp2=config_df.loc[j]["assembly_version"].lower() 
                   
        # "maps/hg38.transcriptID_geneID_map.txt"
        # ENST00000456328	ENSG00000223972	DDX11L1	transcribed_unprocessed_pseudogene
        # ENST00000450305	ENSG00000223972	DDX11L1	transcribed_unprocessed_pseudogene
        # ENST00000488147	ENSG00000227232	WASH7P	unprocessed_pseudogene
        # ENST00000619216	ENSG00000278267	MIR6859-1	miRNA
        # ENST00000473358	ENSG00000243485	MIR1302-2HG	lncRNA

        # "maps/hg38.geneID_geneName_geneType_map.txt"
        # ENSG00000223972	DDX11L1	transcribed_unprocessed_pseudogene
        # ENSG00000227232	WASH7P	unprocessed_pseudogene
        # ENSG00000278267	MIR6859-1	miRNA
        # ENSG00000243485	MIR1302-2HG	lncRNA
        # ENSG00000284332	MIR1302-2	miRNA
        

        sp1_tID_gID_Name_Type = transcript_map_to_dict('maps/%s.transcriptID_geneID_map.txt'%(sp1))
        sp2_tID_gID_Name_Type = transcript_map_to_dict('maps/%s.transcriptID_geneID_map.txt'%(sp2))
        
        for f in ["exons", "genes"]:
            print("\r\t{} {} to {}... finding orthologues".format(i, f, j), end="")
            overlaps = "overlaps/%sto%s.%s.overlap"%(sp1, sp2, f)
            genes = {}
            for line in open(overlaps, 'r'):
                
                # read line overlap between sp1 and sp2
                line = line.strip().split("\t")
                geneA = line[3]
                exonA = ",".join([line[0], line[1], line[2]])
                geneB = line[9]
                exonB = ",".join([line[6], line[7], line[8]]) 

                # calculate overlapping %
                lA = int(line[2])-int(line[1]) #length geneA
                lB = int(line[8])-int(line[7]) #length geneB
                o = int(line[12])              #overlap

                rA = (o*100)/lA                    #overlap of A to B
                rB = 0 if lB == 0 else (o*100)/lB  #overlap of B to A

                # Harmonic mean
                try:
                    hm = 2/((1/rA)+(1/rB))
                except ZeroDivisionError:
                    hm = 0
                
                # Gene_Name and Gene_Type
                geneAtype = "none" if geneA == "." else sp1_tID_gID_Name_Type[geneA]["gene_type"]                
                geneBtype = "none" if geneB == "." else sp2_tID_gID_Name_Type[geneB]["gene_type"]

                # Add sp1 gene overlaps to dictionary
                if not geneA in genes:
                    genes[geneA] = {"overlaps": [],
                                    "genetype": [],
                                    "geneAtype": geneAtype,
                                    "exons": {}}

                if not geneB in genes[geneA]["overlaps"]:
                    genes[geneA]["overlaps"].append(geneB)
                    genes[geneA]["genetype"].append(geneBtype)

                # Add sp1 exon overlaps to dictionary
                if not exonA in genes[geneA]["exons"]:
                    genes[geneA]["exons"][exonA] = {}

                if not geneB in genes[geneA]["exons"][exonA]:
                    genes[geneA]["exons"][exonA][geneB] = [rA, rB, hm]                    

                # Update exonB if new exon have highr A>B overlap
                # NOTE ---> and higher B>A overlap? change % of overlap with HM ??
                if rA > genes[geneA]["exons"][exonA][geneB][0]:
                    genes[geneA]["exons"][exonA][geneB] = [rA, rB, hm]                    
                
                
            # Summarize overlaps into output
            output_file = open('orthology/%sto%s.%s'%(sp1, sp2, f), "w")
            for geneA in genes:
                overlaps = genes[geneA]["overlaps"]
                geneAtype = genes[geneA]["geneAtype"]
                geneBtype = genes[geneA]["genetype"]
                number_of_exons = []
                percentageA = []
                percentageB = []
                for geneB in overlaps:
                    n = 0
                    mA, mB = 0, 0
                    for exon in genes[geneA]["exons"]:
                        if geneB in genes[geneA]["exons"][exon]:
                            n+=1
                            # %overlap A>B (rA)
                            mA+=genes[geneA]["exons"][exon][geneB][0]
                            mB+=genes[geneA]["exons"][exon][geneB][1]
                    number_of_exons.append(str(n))
                    percentageA.append(str(mA/n))
                    percentageA.append(str(mB/n))
                
                percentageA = [ int(float(x)) for x in percentageA ]
                percentageB = [ int(x) for x in percentageB ]
                number_of_exons = [ int(x) for x in number_of_exons ]
    
                percentageA, number_of_exons, overlaps, geneBtype = zip(*sorted(zip(percentageA, number_of_exons, overlaps, geneBtype), reverse=True))
                percentageA = [ str(x) for x in percentageA ]
                number_of_exons = [ str(x) for x in number_of_exons ]
                new_line = "\t".join([geneA, geneAtype, ",".join(overlaps), ",".join(number_of_exons), ",".join(percentageA), ",".join(geneBtype)])
                output_file.write(new_line+"\n")
            output_file.close()
             
            print("\r\t{} {} to {}... done               ".format(i, f, j))

## Classification
##############################Biotypes##############################
# none									none
# stringtie								none
# protein_coding						pc
# 3prime_overlapping_ncRNA				lncRNA
# antisense								lncRNA
# bidirectional_promoter_lncRNA			lncRNA
# lincRNA								lncRNA
# non_coding							lncRNA
# processed_transcript					lncRNA
# sense_intronic						lncRNA
# sense_overlapping						lncRNA
# snoRNA								sncRNA
# snRNA									sncRNA
# sRNA									sncRNA
# misc_RNA								sncRNA
# rRNA									sncRNA
# scaRNA								sncRNA
# TEC 									other
# transcribed_processed_pseudogene 		other
# transcribed_unitary_pseudogene		other
# transcribed_unprocessed_pseudogene	other
# TR_C_gene								other
# TR_D_gene								other
# TR_J_gene								other
# TR_J_pseudogene						other
# TR_V_gene								other
# TR_V_pseudogene						other
# unitary_pseudogene					other
# unprocessed_pseudogene				other
# IG_C_gene								other
# IG_D_gene								other
# IG_J_gene								other
# IG_J_pseudogene						other
# IG_V_gene								other
# IG_V_pseudogene						other
# pseudogene 							other
# ribozyme								other
# polymorphic_pseudogene				other
# processed_pseudogene 					other
####################################################################

check_folder("classification")
print("Classifying orthologues...")

for i in config_df.index:
    for j in config_df.index:
        
        if i == j: continue
    
        print("\r\t{} to {}... classifying".format(i, j), end="")
        sp1=config_df.loc[i]["assembly_version"].lower()
        sp2=config_df.loc[j]["assembly_version"].lower()  
        output_file = open("classification/{}to{}.classification".format(sp1, sp2), 'w')
        
        exons_orthologs = "orthology/{}to{}.exons".format(sp1, sp2)
        genes_orthologs = "orthology/{}to{}.genes".format(sp1, sp2)

        exons = {}
        genes = {}
        types = {}

        for line in open(exons_orthologs, 'r'):
            line = line.strip().split("\t")
            geneM, geneT, ortho, nexon, pcent, btype = parse_orthologs(line)

            counts = count_classes(btype)
            c = classification(counts)
            exons[geneM] = [counts, c, ortho]
            types[geneM] = geneT

        for line in open(genes_orthologs, 'r'):
            line = line.strip().split("\t")
            geneM, geneT, ortho, nexon, pcent, btype = parse_orthologs(line)

            counts = count_classes(btype)
            c = classification(counts)
            genes[geneM] = [counts, c, ortho]
            types[geneM] = geneT

        for gene in exons:
            geneT = types[gene]
            e_counts = exons[gene][0]
            e_class_ = exons[gene][1]
            e_ortho  = exons[gene][2]
            try:
                g_counts = genes[gene][0]
                g_class_ = genes[gene][1]
                g_ortho = genes[gene][2]
            except KeyError:
                g_counts = "."
                g_class_ = "."
                g_ortho = "."

            output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(gene, geneT, e_class_, g_class_, ";".join(e_ortho), e_counts, ";".join(g_ortho), g_counts)+"\n")
        output_file.close()
        print("\r\t{} to {}... done       ".format(i, j))
	
## Merge exon/gene results
## Merge options
# CLASS
#  class1:   eclass == gclass (one to one)
#  class2:   eclass == gclass (many to many)
#  class3:   eclass <= gclass & ! eclass >= glcass
#  class4: ! eclass <= gclass &   eclass >= glcass (This case should not be possible)
#  class5:   eclass == "." & gclass != "." | eclass != "." & gclass == "."
#  class6:   eclass != gclass | eclass == gclass == "."

for i in config_df.index:
    for j in config_df.index:
        if i == j: continue
        print("\r\t{} to {}... assigning class".format(i, j), end="")
        sp1=config_df.loc[i]["assembly_version"].lower()
        sp2=config_df.loc[j]["assembly_version"].lower()
        
        #Read classification data for sp1tosp2
        data = pd.read_csv("./classification/%sto%s.classification"%(sp1, sp2), sep="\t", header=None,
                           names=["gene", "sp1_btype", "eclass", "gclass", "enames", "ecounts", "gnames", "gcounts"]) 

        #Add class and subclass columns to data
        data["class"] = "none"
        data["subclass"] = "none"
        data.set_index("gene", inplace = True)

        #Assign class and subclass
        data=assign_class(data)
        data=assign_subclass(data)
        
        #print(data.loc[data["subclass"]=="subclass1"])
        data.to_csv("./classification/%sto%s.classification.classes"%(sp1, sp2), sep="\t", index=True)
        print("\r\t{} to {}... done            ".format(i, j))
		
### RESULTS
input_path = './classification'
files = [x for x in listdir(input_path) if x.endswith('.classes')]

sps = {}
for input_file in files:
    #Names of sps being analyzed
    org_source = input_file.split('.')[0].split('to')[0]
    org_destin = input_file.split('.')[0].split('to')[1]
    
    #Read orthology results from sp1 to sp2
    df = pd.read_csv(path.join(input_path,input_file),
                     sep='\t',
                     skiprows=1,
                     names=['gene_source', 'sp1_btype', 'e_class','g_class', 'e_names', 'e_counts', 'g_names', 'g_counts', 'class', 'subclass'])
    
    #Count total genes in sp1 and total genes lifted to sp2
    genes_org_source = len(open("BEDs/"+org_source+".genes.bed", 'r').readlines())    
    exons_lifted = set([x.strip().split("\t")[3] for x in  open("liftovers/"+org_source+"to"+org_destin+".exons.liftover", 'r').readlines()])
    genes_lifted = set([x.strip().split("\t")[3] for x in  open("liftovers/"+org_source+"to"+org_destin+".genes.liftover", 'r').readlines()])
    
    #Initialize gene_info dictionary and add total number of genes in sp1
    if not org_source in sps:
        sps[org_source] = {"total_genes": genes_org_source}
    
    #Add information of lifted genes to sp2
    if not org_destin in sps[org_source]:
        
        gene_info = df.set_index('gene_source').T.to_dict()
        sps[org_source][org_destin] = {"total_genes": genes_org_source,
                                       "total_exons_lifted": exons_lifted,
                                       "total_genes_lifted": genes_lifted,
                                       "one_to_one": 0,
                                       "one_to_many": 0,
                                       "many_to_many": 0,
                                       "gene_info": gene_info}
        
gene_type = {}
for input_file in [x for x in os.listdir("./maps") if "geneName" in x]:
    sp = input_file.split(".")[0]
    df_geneType = pd.read_csv(path.join("./maps/",input_file),
                     sep='\t',
                     names=['gene_id', 'gene_name', 'gene_type'])
    gt = df_geneType.set_index('gene_name').T.to_dict()
    gene_type[sp] = gt
	
## Find clusters
# Classification input files    
input_path = './classification'
classification_files = [x for x in listdir(input_path) if x.endswith('.classes')]

# Generate clusters
ALL_GENES = []
EDGES_EXON = []
EDGES_GENE = []

for file in classification_files:

    #Species being analyzed
    org_source = file.split('.')[0].split('to')[0]
    org_destin = file.split('.')[0].split('to')[1]

    #Read pd's with classification
    df = pd.read_csv("./classification/%s"%(file), sep="\t", header=1, names=["gene", "sp1_btype", "eclass", "gclass", "enames", "ecounts", "gnames", "gcounts", "class", "subclass"])
    df = df[['gene','sp1_btype', 'enames','gnames']]

    #Add genes to list
    ALL_GENES += [(org_source, x, gene_type[org_source][x]["gene_type"]) for x in list(df['gene'].unique())]
    ALL_GENES = sorted(set(ALL_GENES))

    #Prepare df's
    df['enames'] = df['enames'].apply(lambda x:[y for y in x.split(';') if y != '.'])
    df['gnames'] = df['gnames'].apply(lambda x:[y for y in x.split(';') if y != '.'])

    #Create edges (gene connections/predictions)
    for row in df.itertuples():

        gene_source = row.gene
        sp1_biotype = row.sp1_btype
        for gene_exon in row.enames:
            sp2_biotype = gene_type[org_destin][gene_exon]["gene_type"]
            edge = ((org_source, gene_source, sp1_biotype), (org_destin, gene_exon, sp2_biotype)) 
            EDGES_EXON.append(edge)

        for gene_gene in row.gnames:
            edge = ((org_source, gene_source, sp1_biotype), (org_destin, gene_gene, sp2_biotype)) 
            EDGES_GENE.append(edge)

# Create network (clusters)
G_exon = nx.DiGraph()
G_exon.add_nodes_from(ALL_GENES)
G_exon.add_edges_from(EDGES_EXON)
G_exon_undirected = G_exon.to_undirected()
components_exon = sorted(list(G_exon_undirected.subgraph(c) for c in nx.connected_components(G_exon_undirected)),key=lambda x:-len(x))
print('EXON graph: %d nodes, %d edges, %d connected_components'%(nx.number_of_nodes(G_exon),
    nx.number_of_edges(G_exon),len(components_exon)),file=sys.stderr)

G_gene = nx.DiGraph()
G_gene.add_nodes_from(ALL_GENES)
G_gene.add_edges_from(EDGES_GENE)
G_gene_undirected = G_gene.to_undirected()
components_gene = sorted(list(G_gene_undirected.subgraph(c) for c in nx.connected_components(G_gene_undirected)),key=lambda x:-len(x))
print('GENE graph: %d nodes, %d edges, %d connected_components'%(nx.number_of_nodes(G_gene),
    nx.number_of_edges(G_gene),len(components_gene)),file=sys.stderr)

exon_to_comp_idx, comp_idx_exon = generate_components_stats(components_exon, G_exon)
gene_to_comp_idx, comp_idx_gene = generate_components_stats(components_gene, G_gene)

## Classify clusters
#For all species
get_clusters_stats(components_exon, G_exon, "exon")
get_clusters_stats(components_gene, G_gene, "gene")

#For 1 vs 1 species
species = list(sps.keys())
if len(species) > 2:
	for sp1_idx in range(0,len(species)):
		for sp2_idx in range(0,len(species)):
			if sp1_idx <= sp2_idx: continue
			sp1 = species[sp1_idx]
			sp2 = species[sp2_idx]

			G_exon = nx.DiGraph()
			G_exon.add_nodes_from([node for node in ALL_GENES if node[0] in [sp1, sp2]])
			G_exon.add_edges_from([edge for edge in EDGES_EXON if edge[0][0] in [sp1, sp2] and  edge[1][0] in [sp1, sp2]])
			G_exon_undirected = G_exon.to_undirected()
			components_exon = sorted(list(G_exon_undirected.subgraph(c) for c in nx.connected_components(G_exon_undirected)),key=lambda x:-len(x))

			G_gene = nx.DiGraph()
			G_gene.add_nodes_from([node for node in ALL_GENES if node[0] in [sp1, sp2]])
			G_gene.add_edges_from([edge for edge in EDGES_GENE if edge[0][0] in [sp1, sp2] and  edge[1][0] in [sp1, sp2]])
			G_gene_undirected = G_gene.to_undirected()
			components_gene = sorted(list(G_gene_undirected.subgraph(c) for c in nx.connected_components(G_gene_undirected)),key=lambda x:-len(x))

			exon_to_comp_idx, comp_idx_exon = generate_components_stats(components_exon, G_exon)
			gene_to_comp_idx, comp_idx_gene = generate_components_stats(components_gene, G_gene)

			get_clusters_stats(components_exon, G_exon, "exon", species_names=[sp1, sp2])
			get_clusters_stats(components_gene, G_gene, "gene", species_names=[sp1, sp2])
        
### Plotting results
check_folder("plots")
subprocess.call("Rscript ./scripts/ConnectOR_GenePlots.R", shell=True)