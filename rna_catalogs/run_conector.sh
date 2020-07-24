#!/usr/bin/bash
# By Migdal 07/2021
# This scripts run ConectOR on several annotation catalogues from human, mouse,
# and zebrafish. The reference catalogues are for human GENCODE v34, mouse
# GENCODE release M24 and for zebrafish ENSEMBL release 100.
# For complete list of catalogues see ref/{human,mouse,zebrafish} directories.
# Each of those catalogues is compared with reference catlogue from other organism.

# Download reference catalogues
if [ ! -d ref ]; then
  mkdir -p ref/human ref/mouse ref/zebrafish
  cd ref/human
  wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.annotation.gtf.gz && mv gencode.v34.annotation.gtf.gz GENCODE_v34.gtf.gz
  wget http://mitranscriptome.org/download/mitranscriptome.gtf.tar.gz && tar -xzf mitranscriptome.gtf.tar.gz && mv mitranscriptome.gtf/mitranscriptome.v2.gtf.gz ./MiTranscriptome_v2.gtf.gz && rm -r  mitranscriptome.gtf.tar.gz mitranscriptome.gtf # http://mitranscriptome.org
  # wget http://www.noncode.org/datadownload/NONCODEv5_human_hg38_lncRNA.gtf.gz # http://www.noncode.org
  # wget http://ccb.jhu.edu/chess/data/chess2.2.gtf.gz # http://ccb.jhu.edu/chess/ -- BedToDict(sp1, f, sp1_geneMap) throws KeyError: 'LOC107985730'
  # BIGTranscriptome by courtesy of Basia bigtrans.lncRNAs.hg38.gene.gtf -- no transcript entries only exons, `generate_maps` can consume only gtf's with transcript entries
  # FANTOM CAT by courtesy of Basia fantomCat.lncRNAs.hg38.gene.gtf -- no transcript entries only exons, `generate_maps` can consume only gtf's with transcript entries
  cd ..
  cd mouse
  wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
  # wget http://www.noncode.org/datadownload/NONCODEv5_mouse_mm10_lncRNA.gtf.gz
  cd ..
  cd zebrafish
  # NONCODE annotation prepared by carlos was used new_danrer10.gtf
  wget ftp://ftp.ensembl.org/pub/release-100/gtf/danio_rerio/Danio_rerio.GRCz11.100.chr.gtf.gz
  cd ../..
fi

# Run ConnectOR for two catalogues
function runConnectOR {
  specie1=${1}
  assembly_version_sp1=${2}
  annotation_sp1=${3}
  specie2=${4}
  assembly_version_sp2=${5}
  annotation_sp2=${6}

  # crete work environment
  home=`pwd`
  ref1=`basename ${annotation_sp1%.gz} .gtf`
  ref2=`basename ${annotation_sp2%.gz} .gtf`
  workdir=${ref1}_${ref2}
  mkdir ${workdir} && cd ${workdir}
  ln -s ${home}/dictionaries.json ${home}/ConnectOR.py ${home}/scripts ${home}/ref/${specie1}/${annotation_sp1} ${home}/ref/${specie2}/${annotation_sp2} .

  # write config
  echo -e "specie\tassembly_version\tannotation\tchainmap\n" > config
  echo -e ${specie1}"\t"${assembly_version_sp1}"\t"${annotation_sp1}"\t\n" >> config
  echo -e ${specie2}"\t"${assembly_version_sp2}"\t"${annotation_sp2}"\t" >> config

  # remove all files from /
  python ConnectOR.py
  title="${ref2} ConnectOR predictions"
  subtitle="${specie1} to ${specie2}"
  echo Rscript scripts/plot_classes.R ${workdir}/classification/ "${title}" "${subtitle}" >> ${home}/plot.R
  cd ${home}  
}

# human to mouse
sp1="human"
assembly_version_sp1="hg38"
sp1_ref="GENCODE_v34.gtf.gz"
sp2="mouse"
assembly_version_sp2="mm10"
for sp2_ref in ref/${sp2}/*; do
  sp2_ref=`basename ${sp2_ref}`
  runConnectOR ${sp1} ${assembly_version_sp1} ${sp1_ref} ${sp2} ${assembly_version_sp2} ${sp2_ref}
done

# human to zebrafish
sp1="human"
assembly_version_sp1="hg38"
sp1_ref="GENCODE_v34.gtf.gz"
sp2="zebrafish"
assembly_version_sp2="danrer11"
for sp2_ref in ref/${sp2}/*; do
  sp2_ref=`basename ${sp2_ref}`
  runConnectOR ${sp1} ${assembly_version_sp1} ${sp1_ref} ${sp2} ${assembly_version_sp2} ${sp2_ref}
done

# mouse to human
sp1="mouse"
assembly_version_sp1="mm10"
sp1_ref="GENCODE_vM25.gtf.gz"
sp2="human"
assembly_version_sp2="hg38"
for sp2_ref in ref/${sp2}/*; do
  sp2_ref=`basename ${sp2_ref}`
  runConnectOR ${sp1} ${assembly_version_sp1} ${sp1_ref} ${sp2} ${assembly_version_sp2} ${sp2_ref}
done

# mouse to zebrafish
sp1="mouse"
assembly_version_sp1="mm10"
sp1_ref="GENCODE_vM25.gtf.gz"
sp2="zebrafish"
assembly_version_sp2="danrer11"
for sp2_ref in ref/${sp2}/*; do
  sp2_ref=`basename ${sp2_ref}`
  runConnectOR ${sp1} ${assembly_version_sp1} ${sp1_ref} ${sp2} ${assembly_version_sp2} ${sp2_ref}
done

# zebrafish to human
sp1="zebrafish"
assembly_version_sp1="danrer11"
sp1_ref="Ensembl_r100.gtf.gz"
sp2="human"
assembly_version_sp2="hg38"
for sp2_ref in ref/${sp2}/*; do
  sp2_ref=`basename ${sp2_ref}`
  runConnectOR ${sp1} ${assembly_version_sp1} ${sp1_ref} ${sp2} ${assembly_version_sp2} ${sp2_ref}
done

# zebrafish to mouse
sp1="zebrafish"
assembly_version_sp1="danrer11"
sp1_ref="Ensembl_r100.gtf.gz"
sp2="mouse"
assembly_version_sp2="mm10"
for sp2_ref in ref/${sp2}/*; do
  sp2_ref=`basename ${sp2_ref}`
  runConnectOR ${sp1} ${assembly_version_sp1} ${sp1_ref} ${sp2} ${assembly_version_sp2} ${sp2_ref}
done
