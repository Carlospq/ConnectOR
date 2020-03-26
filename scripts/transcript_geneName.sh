#!/bin/bash
# Asseign a gene name to each transcript from the novel assembly
# based on the exonic overlap.
# Non overlapping trancripts will be flagged as "NOVEL"

project=$1
bedtools_version=$2
reference_assembly=$3


module add UHTS/Analysis/BEDTools/$bedtools_version
out_dir=$project/sleuth/geneNameMap
mkdir -p $out_dir

#extract transcript_id - gene_name from reference assembly
echo "target_id gene_name" | awk '{OFS="\t"; print $1,$2}' > $out_dir/transcriptID_geneName_map.txt
awk '$3=="transcript"' $reference_assembly |\
 perl -pe 's|.*transcript_id\s"(.*?)";.*gene_name\s"(.*?)";.*|\1 \2|' |\
 sed 's/ /\t/g' >> $out_dir/transcriptID_geneName_map.txt
