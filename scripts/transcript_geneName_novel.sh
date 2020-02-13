#!/bin/bash
# Asseign a gene name to each transcript from the novel assembly
# based on the exonic overlap.
# Non overlapping trancripts will be flagged as "NOVEL"
#
# v1.2(16/09/2019): Extra step added to assigne geneName to transcripts with less than 50% overlap

project=$1
bedtools_version=$2
novel_assembly=$3
reference_assembly=$4


module add UHTS/Analysis/BEDTools/$bedtools_version
out_dir=$project/sleuth/geneNameMap
mkdir -p $out_dir

#extract transcript_id - gene_name from reference assembly
awk '$3=="transcript"' $reference_assembly |\
 perl -pe 's|.*transcript_id\s"(.*?)";.*gene_name\s"(.*?)";.*|\1 \2|' > $out_dir/tr_gn_map_1.txt

#exctract exons from novel assembly without gene_name; output in BED format
grep -v gene_name $novel_assembly |\
 awk '$3=="exon"' |\
 perl -pe 's|(chr.*?)\s.*?(\d+)\s(\d+).*?([+-]).*transcript_id\s"(.*?)";.*?exon_number\s.*?(\d+).*?;.*|\1 \2 \3 \5 \6 \4|' |\
 sed 's/ /\t/g' > $out_dir/exons_wo_geneName.bed

#extract exons from reference assembly (split protein_coding and others if "gene_type" found in reference assembly)
genetype=0
if [[ $(head -n 100000 $reference_assembly | grep transcript_type) ]]; then
 genetype=1;
fi

if [[ $genetype == 1 ]]; then

 awk '$3=="exon"' $reference_assembly |\
 grep "protein_coding" | grep "chr" |\
 perl -pe 's|(chr.*?)\s.*?(\d+)\s(\d+).*?([+-]).*gene_name\s"(.*?)";.*?exon_number\s.*?(\d+).*?;.*|\1 \2 \3 \5 \6 \4|' |\
 sed 's/ /\t/g' > $out_dir/reference_assembly_protein_exons.bed 

 awk '$3=="exon"' $reference_assembly |\
 grep -v "protein_coding" | grep "chr" |\
 perl -pe 's|(chr.*?)\s.*?(\d+)\s(\d+).*?([+-]).*gene_name\s"(.*?)";.*?exon_number\s.*?(\d+).*?;.*|\1 \2 \3 \5 \6 \4|' |\
 sed 's/ /\t/g' > $out_dir/reference_assembly_lnc_exons.bed

else

 awk '$3=="exon"' $reference_assembly | grep chr |\
 perl -pe 's|(chr.*?)\s.*?(\d+)\s(\d+).*?([+-]).*gene_name\s"(.*?)";.*?exon_number\s.*?(\d+).*?;.*|\1 \2 \3 \5 \6 \4|' |\
 sed 's/ /\t/g' > $out_dir/reference_assembly_exons.bed

fi

#intersectBed between SINERGIA & GENCODE exons (first with protein_exons, then with lnc_exons)
if [[ $genetype == 1 ]]; then
 intersectBed -s -loj -f 0.5 -a $out_dir/exons_wo_geneName.bed -b $out_dir/reference_assembly_protein_exons.bed > $out_dir/intersection_pc.txt
 intersectBed -s -loj -a $out_dir/exons_wo_geneName.bed -b $out_dir/reference_assembly_protein_exons.bed > $out_dir/intersection_pc_1nt.txt
 python scripts/filter_lnc_overlaps.py $out_dir/intersection_pc.txt > $out_dir/exons_wo_geneName_lnc.bed
 python scripts/filter_lnc_overlaps.py $out_dir/intersection_pc_1nt.txt > $out_dir/exons_wo_geneName_lnc_1nt.bed
 intersectBed -s -loj -f 0.5 -a $out_dir/exons_wo_geneName_lnc.bed -b $out_dir/reference_assembly_lnc_exons.bed > $out_dir/intersection_lnc.txt
 intersectBed -s -loj -a $out_dir/exons_wo_geneName_lnc_1nt.bed -b $out_dir/reference_assembly_lnc_exons.bed > $out_dir/intersection_lnc_1nt.txt
else
 intersectBed -s -loj -f 0.5 -a $out_dir/exons_wo_geneName.bed -b $out_dir/reference_assembly_exons.bed > $out_dir/intersection.txt
 intersectBed -s -loj -a $out_dir/exons_wo_geneName.bed -b $out_dir/reference_assembly_exons.bed > $out_dir/intersection_1nt.txt
fi

#generate transcript - gene_name map from intersection (for .5% overlap and 1nt overlap)
# .05%
if [[ $genetype == 1 ]]; then
 awk '$10!="."' $out_dir/intersection_pc.txt > $out_dir/intersection_pc_2.txt
 python scripts/mapping.py $out_dir/intersection_pc_2.txt > $out_dir/tr_gn_map_2.txt
 python scripts/mapping.py $out_dir/intersection_lnc.txt >> $out_dir/tr_gn_map_2.txt
else
 python scripts/mapping.py $out_dir/intersection.txt > $out_dir/tr_gn_map_2.txt
fi

# 1nt
if [[ $genetype == 1 ]]; then
 awk '$10!="."' $out_dir/intersection_pc_1nt.txt > $out_dir/intersection_pc_2_1nt.txt
 python scripts/mapping.py $out_dir/intersection_pc_2_1nt.txt > $out_dir/tr_gn_map_2_1nt.txt
 python scripts/mapping.py $out_dir/intersection_lnc_1nt.txt >> $out_dir/tr_gn_map_2_1nt.txt
else
 python scripts/mapping.py $out_dir/intersection_1nt.txt > $out_dir/tr_gn_map_2_1nt.txt
fi

#change non overlapping .05 transcriptID-geneNames with 1nt transcriptID-geneNames overlapping
python scripts/change_novel_geneName.py $out_dir/tr_gn_map_2.txt $out_dir/tr_gn_map_2_1nt.txt > $out_dir/tr_gn_map_2_1nt_2.txt


#merge all mapps in one file
echo "target_id gene_name" | awk '{OFS="\t"; print $1,$2}' > $out_dir/transcriptID_geneName_map.txt
cat $out_dir/tr_gn_map_1.txt $out_dir/tr_gn_map_2.txt >> $out_dir/transcriptID_geneName_map.txt

echo "target_id gene_name" | awk '{OFS="\t"; print $1,$2}' > $out_dir/transcriptID_geneName_map_1nt.txt
cat $out_dir/tr_gn_map_1.txt $out_dir/tr_gn_map_2_1nt_2.txt >> $out_dir/transcriptID_geneName_map_1nt.txt



