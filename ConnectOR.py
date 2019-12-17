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


# In[34]:


### VARIABLES
# minMatch liftOver required to mapp to new region
try:
    minMatch=sys.argv[1]
except IndexError:
    minMatch=50

print(minMatch)
# In[2]:


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
    df = dict(pd.read_csv(file_name, delimiter = dl))
    return(df)

def parse_orthologs(line):
    geneM = line[0]
    ortho = line[1].split(",")
    nexon = line[2].split(",")
    pcent = line[3].split(",")
    btype = line[4].split(",")

    return(geneM, ortho, nexon, pcent, btype)

def count_classes(btype):
    tmp = []
    for type_ in btype:
        if type_ in none:
            tmp.append("none")
        elif type_ in stringtie:
            tmp.append("stringtie")
        elif type_ in lncRNA:
            tmp.append("lncRNA")
        elif type_ in sncRNA:
            tmp.append("sncRNA")
        elif type_ in pc:
            tmp.append("pc")
        else:
            tmp.append("other")
    class_ = [tmp.count("none"), tmp.count("lncRNA"), tmp.count("sncRNA"), tmp.count("pc"), tmp.count("other"), tmp.count("stringtie")]

    return(class_)

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
        c = "lncRNA_PC"
    if classes[0]>=0 and classes[1]>=1 and classes[2]==0 and classes[3]==1 and classes[4]>=1 and classes[5]==0:
        c = "lncRNA_other"
    if classes[0]>=0 and classes[1]>=1 and classes[2]==0 and classes[3]>=1 and classes[4]==0 and classes[5]==0:
        c = "lncRNA_PC"

    return(c)

def change_default(feature, value):
    default_values = config_df.loc[i]["default"].split("|")
    if feature == "gtf": 
        n = 0
    else:
        n = 1
    default_values[n] = value
    return("|".join(default_values))
    


# In[3]:


### Dictionaris
with open('dictionaries.json') as f:
  dictionaries = json.load(f)

#print("dictionaries: ", list(dictionaries.keys()), "\n")
#for k in dictionaries:
#    print(k, list(dictionaries[k].keys()))
    
#print(dictionaries["chain_maps"]["danrer10"])


# ### Read config file

# In[4]:


config_df = read_config("./config")
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

# In[5]:


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

# In[6]:


print('Downloading GTFs')
for i in config_df.index:
    if config_df.loc[i]["default"].split("|")[0] == "False": 
        print("\tNo default GTF needed for %s.. Skipping"%(i))
        continue    
    sp_v = config_df.loc[i]["assembly_version"].lower()
    download_default_files(dictionaries["gtfs_ensembl_r98"][sp_v], "GTFs")        


# ### Download chainmaps

# In[7]:


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

# In[19]:


check_folder("maps")
print("Generating transcriptID-geneID & geneID-geneName-geneType maps...")
for i in config_df.index:
    gtf_file = config_df.loc[i]["annotation"]
    sp_v = config_df.loc[i]["assembly_version"].lower()
    if config_df.loc[i]["default"].split("|")[0] == "True":
        gtf_file = "GTFs/"+gtf_file
    print("\t"+gtf_file)
    generate_maps(gtf_file, sp_v)


# In[17]:


print("Transcripts map:")
with open("maps/hg38.transcriptID_geneID_map.txt") as transcripts:
    head = [next(transcripts) for x in range(5)]
for l in head:
    print(l.strip())

print("\nGenes map:")
with open("maps/hg38.geneID_geneName_geneType_map.txt") as genes:
    head = [next(genes) for x in range(5)]
for l in head:
    print(l.strip())


# ### Generate BED with genes/exons from GTF files

# In[24]:


check_folder("BEDs")
print("Generating BED files for exons and genes...")
for i in config_df.index:
    gtf_file = config_df.loc[i]["annotation"]
    sp_v = config_df.loc[i]["assembly_version"].lower()
    if config_df.loc[i]["default"].split("|")[0] == "Treu":
        gtf_file = "GTFs/"+gtf_file
    print("\t"+gtf_file)
    generate_beds(gtf_file, sp_v)


# In[25]:


print("Sorting BED files...")
for i in config_df.index:
    sp_v = config_df.loc[i]["assembly_version"].lower()
    bed_sort(sp_v)


# ### LifOver exons/genes

# In[29]:


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
            os.system("./liftOver -minMatch=0.{} {} {} {} {}".format(minMatch, oldFile, map_chain, newFile, unMapped))
            print("\r\t{} {} to {}... done   ".format(i, feature, j))


# ### Intersect LiftOvers

# In[30]:


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


# ### Parse Overlaps

# In[31]:


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
        geneID_Name_Type = gene_map_to_dict('maps/%s.geneID_geneName_geneType_map.txt'%(sp2))

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
                lA = int(line[2])-int(line[1])
                lB = int(line[8])-int(line[7])
                o = int(line[12])

                rA = (o*100)/lA
                rB = 0 if lB == 0 else (o*100)/lB

                # Harmonic mean
                try:
                    hm = 2/((1/rA)+(1/rB))
                except ZeroDivisionError:
                    hm = 0

                # Gene_Name and Gene_Type from sp2 gene 
                gene_name = "." if geneB == "." else geneID_Name_Type[geneB]["gene_name"] #not needed so far
                geneBtype = "none" if geneB == "." else geneID_Name_Type[geneB]["gene_type"]
            
                # Add sp1 gene overlaps to dictionary
                if not geneA in genes:
                    genes[geneA] = {"overlaps": [],
                                    "genetype": [],
                                    "exons": {}}

                if not geneB in genes[geneA]["overlaps"]:
                    genes[geneA]["overlaps"].append(geneB)
                    genes[geneA]["genetype"].append(geneBtype)

                # Add sp1 exon overlaps to dictionary
                if not exonA in genes[geneA]["exons"]:
                    genes[geneA]["exons"][exonA] = {}

                if not geneB in genes[geneA]["exons"][exonA]:
                    genes[geneA]["exons"][exonA][geneB] = [rA, rB, hm]

                # Exon with highr A>b and B>A overlap 
                # change % of overlap with HM ??
                if rA > genes[geneA]["exons"][exonA][geneB][0]:
                    genes[geneA]["exons"][exonA][geneB] = [rA, rB, hm]

            # Summarize overlaps into output
            output_file = open('orthology/%sto%s.%s'%(sp1, sp2, f), "w")
            for geneA in genes:
                overlaps = genes[geneA]["overlaps"]
                genetype = genes[geneA]["genetype"]
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
    
                percentageA, number_of_exons, overlaps, genetype = zip(*sorted(zip(percentageA, number_of_exons, overlaps, genetype), reverse=True))
                percentageA = [ str(x) for x in percentageA ]
                number_of_exons = [ str(x) for x in number_of_exons ]
                new_line = "\t".join([geneA, ",".join(overlaps), ",".join(number_of_exons), ",".join(percentageA), ",".join(genetype)])
                output_file.write(new_line+"\n")
            output_file.close()
             
            print("\r\t{} {} to {}... done               ".format(i, f, j))


# ### Classification

# In[32]:


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

        none = ["none"]
        stringtie = ["stringtie"]
        pc = ["protein_coding"]
        lncRNA = ["3prime_overlapping_ncRNA", "antisense", "bidirectional_promoter_lncRNA", "lincRNA", "non_coding", "processed_transcript", "sense_intronic", "sense_overlapping", "lncRNA"]
        sncRNA = ["snoRNA", "snRNA", "sRNA", "misc_RNA", "rRNA", "scaRNA", "miRNA"]

        exons = {}
        genes = {}

        for line in open(exons_orthologs, 'r'):
            line = line.strip().split("\t")
            geneM, ortho, nexon, pcent, btype = parse_orthologs(line)

            counts = count_classes(btype)
            c = classification(counts)
            exons[geneM] = [counts, c, ortho]

        for line in open(genes_orthologs, 'r'):
            line = line.strip().split("\t")
            geneM, ortho, nexon, pcent, btype = parse_orthologs(line)

            counts = count_classes(btype)
            c = classification(counts)
            genes[geneM] = [counts, c, ortho]

        for gene in exons:
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

            output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(gene, e_class_, g_class_, ";".join(e_ortho), e_counts, ";".join(g_ortho), g_counts)+"\n")
        output_file.close()
        print("\r\t{} to {}... done       ".format(i, j))


# ### Plots

# In[33]:


check_folder("plots")
print("Plotting results...")

for i in config_df.index:
    for j in config_df.index:
        print("\r\t{} to {}... plotting".format(i, j), end="")
        if i == j: continue
        sp1 = config_df.loc[i]["assembly_version"].lower()
        sp2 = config_df.loc[j]["assembly_version"].lower() 
        input_file = "classification/{}to{}.classification".format(sp1, sp2)
        path = "./"
        call = 'Rscript plots.R %s %s %s %s'%(input_file, sp1, sp2, path)
        #print(call)
        print("\r\t{} to {}... done       ".format(i, j))
        subprocess.call(call, shell=True, executable='/bin/bash') 
        


# In[ ]:





# In[ ]:




