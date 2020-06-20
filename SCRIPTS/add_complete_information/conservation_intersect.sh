home="../..";

conservation_home="${home}/REFERENCE_DATA/conservation";
riboseq_software="${home}/REFERENCE_DATA/riboseq/software";




#### uncompress the alignments into bedgraph format
# http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phastCons100way/hg38.phastCons100way.bw
#$riboseq_software/bigWigToBedGraph ${conservation_home}/hg38.phastCons100way.bw ${conservation_home}/hg38.phastCons100way.bg;

# http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/hg38.phyloP100way.bw
#$riboseq_software/bigWigToBedGraph ${conservation_home}/hg38.phyloP100way.bw ${conservation_home}/hg38.phyloP100way.bg;

input_filename="$1";
filt="${home}/DATA/filtered";
batch_home="$filt/$input_filename";

conservation_output_home="$batch_home/conservation";



# Here we put the phastCons and the phyloP scores for each feature of the uORF in the same row
for element in start stop exon; do

   echo "merging phastCons and phyloP ${element}";

   awk 'BEGIN{FS=OFS="\t"}{a[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15]=a[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15]"\t"$16;}END{ for (i in a) print  i""a[i] }' $conservation_output_home/phastCons100way/uORFs_phastCons100way_conservation_joined.tsv.${element} $conservation_output_home/phyloP100way/uORFs_phyloP100way_conservation_joined.tsv.${element} | sort -k1,1 -k2,2n > $conservation_output_home/all_conservation_by_${element}.tsv;


done;


# merge multiple exons into a unique entry, same for splitted start and stop codons
unifying="$conservation_output_home/unifying";
mkdir -p $unifying;

for element in start stop exon; do
   echo "combining same feature information in ${element}";

   awk 'BEGIN{FS=OFS="\t"}{  \
         a[$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15]+=$16; \
         b[$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15]+=$17; \
         }END{ for (i in a) print  i"\t"a[i]"\t"b[i] }' $conservation_output_home/all_conservation_by_${element}.tsv | sort -k1,1 -k2,2n > $unifying/all_conservation_unified_${element}.tsv;

done;





echo "Joining start, exon and stop in this order";
awk 'BEGIN{FS=OFS="\t"}{a[$1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12]=a[$1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12]"\t"$13"\t"$14"\t"$15;}END{ for (i in a) print  i""a[i] }' $unifying/all_conservation_unified_start.tsv $unifying/all_conservation_unified_exon.tsv $unifying/all_conservation_unified_stop.tsv | sort -k5,5 -k4,4 > $conservation_output_home/all_conservation_by_uORF.tsv;

#??? removing empty exon and stop codon matching??
# cut -f 15,18 --complement $conservation_output_home/all_conservation_by_uORF.tsv > $conservation_output_home/final_all_conservation_by_uORF.tsv;

echo "Remove the column with exon matching reference";
cut -f 15 --complement $conservation_output_home/all_conservation_by_uORF.tsv > $conservation_output_home/final_all_conservation_by_uORF.tsv;

rm $conservation_output_home/all_conservation_by_uORF.tsv;

# we may modify some columns (ID) to make it easier to compare and if I could join those loci that are the same it would also be great

# 7 nucleotide length ( we must add 3 because the stop codon is not considered and we want to consider it)
# 13-21
# chr   uORFtype strand   uORF_id  gene_id  transcript_id  nuc_len  aa_len   extra cdna_seq prot_seq existingStartCodonIDs   start_phastCons100way start_phyloP100way   exon_phastCons100way exon_phyloP100way stop_phastCons100way stop_phyloP100way
# chrX	upstream_uORF	-	uORF-ENST00000612152-4;chrX:100637026-100637040:-	ENSG00000000003	ENST00000612152	12	4	.	ATGGTGAAGTTATAA	MVKL*	.	0.145	1.275	0.267	6.254	0.011	1.952

#using phyloP total
#awk 'BEGIN {FS=OFS="\t";} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13/3, $14, $15/($7+3), $16, $17/3, $18}' \
#      $filt/final_all_conservation_by_uORF.tsv > $filt/final_all_conservation_by_uORF.average.tsv;

##using phyloP average
#awk 'BEGIN {FS=OFS="\t";} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13/3, $14/3, $15/($7+3), $16/($7+3), $18/3, $19/3}' \
#      $conservation_output_home/final_all_conservation_by_uORF.tsv > $conservation_output_home/final_all_conservation_by_uORF.average.tsv;
awk 'BEGIN {FS=OFS="\t";} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12":"$17, $13/3, $14/3, $15/($7+3), $16/($7+3), $18/3, $19/3}' \
      $conservation_output_home/final_all_conservation_by_uORF.tsv > $conservation_output_home/final_all_conservation_by_uORF.average.tsv;


echo "Compressing previous files";
tar -czvf $conservation_output_home/intermediate_uORFs_conservation_files.tar.gz $conservation_output_home/all_conservation_by_* $conservation_output_home/phastCons100way/ $conservation_output_home/phyloP100way/  $conservation_output_home/unifying/;

echo "Compressed now we remove them";
rm -rf $conservation_output_home/all_conservation_by_* $conservation_output_home/phastCons100way/ $conservation_output_home/phyloP100way/  $conservation_output_home/unifying/;


echo "--- DONE ---";
echo "uORFs with conservation data can be found here:";
echo "$conservation_output_home/final_all_conservation_by_uORF.average.tsv";




