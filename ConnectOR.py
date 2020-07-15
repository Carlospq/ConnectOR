#!/usr/bin/env python
# coding: utf-8

# # ConectOR (Conect Orthologue RNAs)

# In[1]:


### IMPORTS
import re
import subprocess
import os, sys
import fileinput
import pandas as pd
from io import StringIO
import json
import wget
import gzip
#from tqdm import tqdm


# In[34]:


### VARIABLES
# minMatch liftOver required to mapp to new region
try:
    minMatch=sys.argv[1]
except IndexError:
    minMatch=50
    
# sys.argv[] does not work properly in jupyter
minMatch = 50

### FUNCTIONS

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
        print("".join(["\r\t",file_name,"... sorting"]), end = '')
        call = "sort -u -k1,1 -k2,2n -o '%s' '%s'"%(file_name,file_name)    
        subprocess.call(call, shell=True)
        print("".join(["\r\t",file_name,"... sorted "]))
    

def merge_bed(spv):
    #This functions assumes BED files are sorted
    files = ["./BEDs/{}.exons.bed".format(spv)]
    genes = {}
    previous_id = ""
    
    for file_name in files:
        print("\r\t"+file_name+"... merging", end="")
        for line in open(file_name, 'r'):
            line=line.strip().split("\t")

            chrom = line[0]
            start = line[1]
            end = line[2]
            gene_name = line[3]
            score = line[4]
            strand = line[5]
            
            iexon = [chrom, start, end, gene_name, score, strand]
            
            if not gene_name in genes:
                genes[gene_name] = [iexon]
                continue

            jexon = genes[gene_name][-1]

            if int(iexon[1]) >= int(jexon[1]) and int(iexon[1]) <= int(jexon[2]):
                
                if int(iexon[2]) >= int(jexon[2]):
                    jexon[2]=iexon[2]

            if int(iexon[1]) > int(jexon[1]) and int(iexon[1]) > int(jexon[2]):
                if not iexon in genes[gene_name]:
                    genes[gene_name].append(iexon)
        
        output_exons = open("./BEDs/{}.exons.bed".format(spv),"w")
        for gene in genes:
            for exon in genes[gene]:
                exon_bed_line = "\t".join(exon)+"\n"
                output_exons.write(exon_bed_line)
        output_exons.close()
        print("\r\t"+file_name+"... merged ")
        
        
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
    tmp = []
    for type_ in btype:
        if type_ in none:
            tmp.append("none")
        elif type_ in lncRNA:
            tmp.append("lncRNA")
        elif type_ in pc:
            tmp.append("pc")
        elif type_ in pseudogene:
            tmp.append("pseudogene")
        elif type_ in non_lifted:
            tmp.append("non_lifted")
        else:
            tmp.append("other")
    class_ = [tmp.count("none"), tmp.count("lncRNA"), tmp.count("pc"), tmp.count("pseudogene"), tmp.count("other"), tmp.count("non_lifted")]

    return(class_)


def classification_tree(exon_info, gene_info):
    #counts=["none", "lncRNA", "pc", "pseudogene", "other", "non_lifted"]
    #var0=atype, var1=btype, var2=counts, var3=ortho  
    exon_dict=dict(("var{0}".format(i),x) for i,x in enumerate(exon_info))
    gene_dict=dict(("var{0}".format(i),x) for i,x in enumerate(gene_info))

    e_counts = exon_dict["var2"]
    g_counts = gene_dict["var2"]
    if exon_dict["var1"][0] == "non_lifted" or gene_dict["var1"][0] == "non_lifted":
        prediction = "non_lifted"
    elif e_counts[2] >= 1:
        prediction = "pc(h.c.)"
    elif e_counts[1] >= 1:
        prediction = "lncRNA(h.c.)"
    elif g_counts[2] >= 1:
        prediction = "pc(l.c.)"
    elif g_counts[2] >= 1:
        prediction = "lncRNA(l.c.)"
    elif g_counts[3] >= 1:
        prediction = "pseudogene(l.c.)"
    elif g_counts[4] >= 1:
        prediction = "other"
    elif g_counts[0] >= 1:
        prediction = "none"

    return(prediction)


def classification(classes):
    c = "other"
    #unique cases
    #  none 			 lncRNA 		   sncRNA 			 pc 			   other 			 stringtie
    if classes[0]==1 and classes[1]==0 and classes[2]==0 and classes[3]==0 and classes[4]==0 and classes[5]==0:
        c = "none"
    if classes[0]>=0 and classes[1]==1 and classes[2]==0 and classes[3]==0 and classes[4]==0 and classes[5]==0:
        c = "lncRNA"
    if classes[0]>=0 and classes[1]==0 and classes[2]==1 and classes[3]==0 and classes[4]==0 and classes[5]==0:
        c = "sncRNA"
    if classes[0]>=0 and classes[1]==0 and classes[2]==0 and classes[3]==1 and classes[4]==0 and classes[5]==0:
        c = "pc"
    if classes[0]>=0 and classes[1]==0 and classes[2]==0 and classes[3]==0 and classes[4]==1 and classes[5]==0:
        c = "other"
    if classes[0]>=0 and classes[1]==0 and classes[2]==0 and classes[3]==0 and classes[4]==0 and classes[5]>0:
        c = "stringtie"

    #multiple cases
    #  none 			 lncRNA 		   sncRNA 			 pc 			  others			stringtie
    if classes[0]>1 and classes[1]==0 and classes[2]==0 and classes[3]==0 and classes[4]==0 and classes[5]==0:
        c = "nones"
    if classes[0]>=0 and classes[1]>1 and classes[2]==0 and classes[3]==0 and classes[4]==0 and classes[5]==0:
        c = "lncRNAs"
    if classes[0]>=0 and classes[1]==0 and classes[2]>1 and classes[3]==0 and classes[4]==0 and classes[5]==0:
        c = "sncRNAs"
    if classes[0]>=0 and classes[1]==0 and classes[2]==0 and classes[3]>1 and classes[4]==0 and classes[5]==0:
        c = "pcs"
    if classes[0]>=0 and classes[1]==0 and classes[2]==0 and classes[3]==0 and classes[4]>1:
        c = "others"

    #dual cases
    #  none 			 lncRNA 		   sncRNA 			 pc 			   other
    if classes[0]>=0 and classes[1]>=1 and classes[2]>=1 and classes[3]==0 and classes[4]==0 and classes[5]==0:
        c = "lncRNA_sncRNA"
    if classes[0]>=0 and classes[1]>=1 and classes[2]==0 and classes[3]==0 and classes[4]==0 and classes[5]>=1:
        c = "lncRNA_stringtie"
    if classes[0]>=0 and classes[1]>=1 and classes[2]==0 and classes[3]>=1 and classes[4]==0 and classes[5]==0:
        c = "lncRNA_pc"
    if classes[0]>=0 and classes[1]>=1 and classes[2]==0 and classes[3]==1 and classes[4]>=1 and classes[5]==0:
        c = "lncRNA_other"
    if classes[0]>=0 and classes[1]>=1 and classes[2]==0 and classes[3]>=1 and classes[4]==0 and classes[5]==0:
        c = "lncRNA_pc"

    return(c)


def change_default(feature, value):
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
                print(set1, set2)
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


### Dictionaris
with open('dictionaries.json') as f:
  dictionaries = json.load(f)

# print("dictionaries: ", list(dictionaries.keys()), "\n")
# for k in dictionaries:
#     print(k, list(dictionaries[k].keys()))
    
# print(dictionaries["chain_maps"]["danrer10"])

### Read config file

config_df = read_config("./config")
#config_df = read_config("./config_local")
config_df.set_index("specie", inplace = True)
config_df["default"] = "False|False"
for i in config_df.index:
    if not config_df.loc[i]["annotation"]:
        default_gtf = dictionaries["gtfs_ensembl_r98"][config_df.loc[i]["assembly_version"].lower()].split("/")[-1]
        config_df.at[i, 'default'] = change_default("gtf", "True")
        config_df.at[i, 'annotation'] = default_gtf
    if not config_df.loc[i]["chainmap"]:
        chainmaps = []
        config_df.at[i, 'default'] = change_default("chainmap", "True")
        for j in config_df.index:
            if i!=j:        
                default_chainmap_path = dictionaries["chain_maps"][config_df.loc[i]["assembly_version"].lower()][config_df.loc[j]["assembly_version"].lower()]
                default_chainmap_name = "chainmaps/"+default_chainmap_path.split("/")[-1]
                chainmaps.append(default_chainmap_name)
        config_df.at[i, 'chainmap'] = ",".join(chainmaps)
config_df


# ### Check config file

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


# ### Download GTF files (if no gtf provided in config)

print('Downloading GTFs')
for i in config_df.index:
    if config_df.loc[i]["default"].split("|")[0] == "False": 
        print("\tNo default GTF needed for %s.. Skipping"%(i))
        continue    
    sp_v = config_df.loc[i]["assembly_version"].lower()
    download_default_files(dictionaries["gtfs_ensembl_r98"][sp_v], "GTFs")            


# ### Download chainmaps

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


# ### Generate maps (transcriptID-geneID & geneID-geneName-geneType) 

check_folder("maps")
print("Generating transcriptID-geneID & geneID-geneName-geneType maps...")
for i in config_df.index:
    print("\r\tGenerating map for "+i+"...", end="")
    gtf_file = config_df.loc[i]["annotation"]
    sp_v = config_df.loc[i]["assembly_version"].lower()
    if config_df.loc[i]["default"].split("|")[0] == "True":
        gtf_file = "GTFs/"+gtf_file
    generate_maps(gtf_file, sp_v)
    print("\r\tGenerating map for "+i+"... done")


# ### Generate BED with genes/exons from GTF files

check_folder("BEDs")
print("Generating BED files for exons and genes...")
for i in config_df.index:
    gtf_file = config_df.loc[i]["annotation"]
    sp_v = config_df.loc[i]["assembly_version"].lower()
    if config_df.loc[i]["default"].split("|")[0] == "True":
        gtf_file = "GTFs/"+gtf_file
    print("\r\tGenerating BEDs for "+i+"...", end="")
    generate_beds(gtf_file, sp_v)
    print("\r\tGenerating BEDs for "+i+"... done")
    
print("Sorting BED files...")
for i in config_df.index:
    sp_v = config_df.loc[i]["assembly_version"].lower()
    bed_sort(sp_v)
    
print("Merging BED files...")
for i in config_df.index:
    sp_v = config_df.loc[i]["assembly_version"].lower()
    merge_bed(sp_v)

# ### LifOver exons/genes

check_folder("liftovers")
features = ["exons", "genes"]
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
            print("\r\t{} {} to {}... mapping".format(i, feature, j), end="")
            oldFile = "BEDs/{}.{}.bed".format(sp_vi, feature)
            newFile = "liftovers/{}to{}.{}.liftover".format(sp_vi, sp_vj, feature)
            unMapped= "liftovers/{}to{}.{}.unmapped".format(sp_vi, sp_vj, feature)
            #print("./liftOver {} {} {} {}".format(oldFile, map_chain, newFile, unMapped))
            os.system("./scripts/liftOver -multiple -minMatch=0.{} {} {} {} {}".format(minMatch, oldFile, map_chain, newFile, unMapped))
            print("\r\t{} {} to {}... done   ".format(i, feature, j))

# ### Intersect LiftOvers

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
            call = '/home/mmigdal/miniconda3/bin/intersectBed -wao -s -a %s -b %s > %s'%(lifover_input, bed_input, output)
            subprocess.call(call, shell=True, executable='/bin/bash')
            print("\r\t{} {} to {}... done         ".format(i, f, j))


# ### Parse Overlaps

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
        
        genes = {}
        sp1_tID_gID_Name_Type = transcript_map_to_dict('maps/%s.transcriptID_geneID_map.txt'%(sp1))
        sp2_tID_gID_Name_Type = transcript_map_to_dict('maps/%s.transcriptID_geneID_map.txt'%(sp2))
        
        for f in ["exons", "genes"]:
            print("\r\t{} {} to {}... finding orthologues".format(i, f, j), end="")
            overlaps = "overlaps/%sto%s.%s.overlap"%(sp1, sp2, f)
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

### Add Unmapped exons/genes

def BedToDict(sp1, feature, geneMap):
    inFile = open("./BEDs/%s.%s.bed"%(sp1, feature), 'r')
    d = {}
    for line in inFile:
        line = line.strip().split("\t")
        GeneName = line[3]
        chrom = line[0]
        start = line[1]
        end = line[2]
        strand = line[5]
        gType = geneMap[GeneName]["gene_type"]
        
        if not GeneName in d:
            d[GeneName] = [chrom, start, end, GeneName, "0", strand]
    inFile.close()
    return(d)


def OrthologyToList(sp1, sp2, feature):
    inFile = open("./orthology/%sto%s.%s"%(sp1, sp2, feature), 'r')
    l = []
    for line in inFile:
        line = line.strip().split("\t")
        GeneName = line[0]
        l.append(GeneName)
    inFile.close()
    return(l)

print("Adding non-lifted genes...")
for i in config_df.index:
    for j in config_df.index:
        
        if i == j: continue
        print("\r\tAdding non-lifted {} genes to results...".format(i, f, j), end="")
        sp1=config_df.loc[i]["assembly_version"].lower()
        sp2=config_df.loc[j]["assembly_version"].lower()
        sp1_geneMap = gene_map_to_dict('maps/%s.geneID_geneName_geneType_map.txt'%(sp1))
        
        for f in ["exons", "genes"]:
            d = BedToDict(sp1, f, sp1_geneMap)
            l = OrthologyToList(sp1, sp2, f)
            #print("before: ", len(d))
            [d.pop(x) for x in l if x in d]
            #print("after: ", len(d))
            
            fileName = "./orthology/"+sp1+"to"+sp2+"."+f
            outFile = open(fileName, 'a')
            for k in d:
                info = d[k]
                new_line = "\t".join([k, info[3], "-", "0", "0", "non_lifted\n"])
                outFile.write(new_line)
            outFile.close()
        print("\r\tAdding non-lifted {} genes to results... done".format(i, f, j))

# ### Classification

for i in config_df.index:
    for j in config_df.index:
        
        if i == j: continue
    
        sp1=config_df.loc[i]["assembly_version"].lower()
        sp2=config_df.loc[j]["assembly_version"].lower()  
        
        exons_orthologs = "orthology/{}to{}.exons".format(sp1, sp2)
        genes_orthologs = "orthology/{}to{}.genes".format(sp1, sp2)
        h = []
        for line in open(exons_orthologs, 'r'):
            line=line.strip().split("\t")
            bt = line[-1].split(",")
            for e in bt:
                if not e in h:
                    h.append(e)
    
for i in set(h):
    print(i)
    
##############################Biotypes##############################
# non_lifted                          non_lifted
# none                                none
# lncRNA                              lncRNA
# protein_coding                      pc
# pseudogene                          pseudogene
# unitary_pseudogene                  pseudogene
# translated_processed_pseudogene     pseudogene
# unprocessed_pseudogene              pseudogene
# transcribed_unprocessed_pseudogene  pseudogene
# transcribed_processed_pseudogene    pseudogene
# transcribed_unitary_pseudogene      pseudogene
# polymorphic_pseudogene              pseudogene
# processed_pseudogene                pseudogene
####################################################################

check_folder("classification")
print("Classifying orthologues...")
    
for i in config_df.index:
    for j in config_df.index:
        
        if i == j: continue
    
        print("\r\t{} to {}... classifying".format(i, j))
        sp1=config_df.loc[i]["assembly_version"].lower()
        sp2=config_df.loc[j]["assembly_version"].lower()  
        output_file = open("classification/{}to{}.classification".format(sp1, sp2), 'w')
        
        exons_orthologs = "orthology/{}to{}.exons".format(sp1, sp2)
        genes_orthologs = "orthology/{}to{}.genes".format(sp1, sp2)

        non_lifted = ["non_lifted"]
        none = ["none"]
        lncRNA = ["lncRNA"]
        pc = ["protein_coding"]
        pseudogene = ["pseudogene", "unitary_pseudogene", "translated_processed_pseudogene", "unprocessed_pseudogene", "transcribed_unprocessed_pseudogene", "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", "polymorphic_pseudogene", "processed_pseudogene"]

        exons = {}
        genes = {}

        for line in open(exons_orthologs, 'r'):
            line = line.strip().split("\t")
            geneM, atype, ortho, nexon, pcent, btype = parse_orthologs(line)
            
            exon_counts = count_classes(btype)
            exons[geneM] = [exon_counts, ortho, atype, btype]
            
        for line in open(genes_orthologs, 'r'):
            line = line.strip().split("\t")
            geneM, atype, ortho, nexon, pcent, btype = parse_orthologs(line)

            gene_counts = count_classes(btype)
            genes[geneM] = [gene_counts, ortho, atype, btype]
      
        for gene in exons:
            e_atype = exons[gene][2]
            e_btype = exons[gene][3]
            e_counts = exons[gene][0]
            e_ortho  = exons[gene][1]
            exon_info = [e_atype, e_btype, e_counts, e_ortho]
            try:
                g_atype = genes[gene][2]
                g_btype = genes[gene][3]
                g_counts = genes[gene][0]
                g_ortho = genes[gene][1]
            except KeyError:
                g_counts = "."
                g_class_ = "."
                g_ortho = "."
            gene_info = [g_atype, g_btype, g_counts, g_ortho]             
            prediction = classification_tree(exon_info, gene_info)

            g_atype = "pseudogene" if "pseudogene" in g_atype else g_atype
            output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(gene, g_atype, prediction, ";".join(e_ortho), ";".join(e_btype), e_counts, ";".join(g_ortho), ";".join(g_btype), g_counts)+"\n")
        output_file.close()
        print("\r\t{} to {}... done       ".format(i, j))
        

