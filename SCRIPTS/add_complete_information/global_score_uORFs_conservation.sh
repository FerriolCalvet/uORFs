
# home="/home/fcalvet/Desktop/uORFs";
home=".";

scripts="$home/SCRIPTS/global_scripts";

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



# we only perform this step if create is equal to 1, otherwise we assume that the input files have already been created
# as we execute this file from the processing_add_features one, this has already been created
create="0";

if [ "$create" == "1" ]; then
   ### START prepare the uORF candidates file for merging with any RP data or conservation using BedTools

   echo "Creation begins";
   for element in start stop exon; do
      echo "Starting ${element}";
   #   element="start";
   #   element="stop";
   #   element="exon";

      uORF_candidates="$batch_home/nh_${element}_merged_$input_filename.tsv";
      uORF_input="$batch_home/nh_m_${element}_input_$input_filename.tsv";

      # rearrange the uORFs file. Use this for elongation. (intersect the whole uORF)
      awk 'BEGIN {FS=OFS="\t";} {print $1, $4-1, $5, $2, $3, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15}' \
            $uORF_candidates > $uORF_input.tmp;

      sort -k1,1 -k2,2n  $uORF_input.tmp > $uORF_input;
      
      rm $uORF_input.tmp;

   done;
   echo "Creation DONE!";
   echo "";

   ### END prepare the uORF candidates file for merging with any RP data or conservation using BedTools
fi





# all the files created in the process of adding conservation features to the uORFs will be stored in this directory
conservation_output_home="$batch_home/conservation";
mkdir -p $conservation_output_home;


# this loop allows us to process the RP coverage data and merge it with the uORFs candidates
for specific in phastCons100way phyloP100way; do
# for specific in phyloP100way; do

#   specific=phyloP100way;
#   specific=phastCons100way;

   echo "Intersection Start ${specific}";

   specificConservationOutput="$conservation_output_home/${specific}";
   mkdir -p $specificConservationOutput;

#   specificRiboSeqInput="${conservation_home}";

   global_accumulated="$conservation_home/hg38.${specific}.bg";

   # this means we are decompressing the phyloP or phastCons file every time we run this script
   echo "Decompressing ${specific}";
   $riboseq_software/bigWigToBedGraph ${conservation_home}/hg38.${specific}.bw ${global_accumulated};
   # rm ${conservation_home}/hg38.${specific}.bw;

   for element in start stop exon; do
   #for element in exon; do
      echo "Starting ${specific} ${element} ";
      #   element="start";
      #   element="stop";
      #   element="exon";

      uORF_input="$batch_home/nh_m_${element}_input_$input_filename.tsv";


      ### START intersect

      # columns need to be adjusted to the specific cases we are handling, this is for initiating
      # intersection between the bedgraph with all the coverage of all files, and the uORF candidates
      # I flipped a and b to have a proper combination and output with the "wao" intersection
      # # bedtools intersect -wao -a $global_accumulated -b $specificRiboSeq/uORFs_input_${specific}.tsv > $specificRiboSeq/uORFs_${specific}_RPcoverage_pre.tsv
#      bedtools intersect -sorted -wao -a $uORF_input -b $global_accumulated > $specificConservationOutput/uORFs_${specific}_conservation_pre.tsv.${element};
#      echo "Bedtools done! ";
#done;

      echo "Starting bedtools intersect";
#      bedtools intersect -wao -a $uORF_input -b $global_accumulated > $specificConservationOutput/uORFs_${specific}_conservation_pre.tsv.${element};

# Apparently, according to the documentation of bedtools intersect, sorted option makes the program faster but also uses more memory
      bedtools intersect -sorted -wao -a $uORF_input -b $global_accumulated > $specificConservationOutput/uORFs_${specific}_conservation_pre.tsv.${element};
      echo "Finishing bedtools intersect";



      awk '{ print $0"\t"$19 * $20 }' $specificConservationOutput/uORFs_${specific}_conservation_pre.tsv.${element} > $specificConservationOutput/uORFs_${specific}_conservation.tsv.${element};
      rm $specificConservationOutput/uORFs_${specific}_conservation_pre.tsv.${element};

      awk 'BEGIN{FS=OFS="\t"}{a[$1"\t"$2+1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15]+=$21;}END{ for (i in a) print  i"\t"a[i] }' $specificConservationOutput/uORFs_${specific}_conservation.tsv.${element} | sort -k1,1 -k2,2n > $specificConservationOutput/uORFs_${specific}_conservation_joined.tsv.${element};

      # process the file to recompute the total coverage of each region, we do it by multiplying the
      # coverage value of an entry times the position with overlap

      rm $specificConservationOutput/uORFs_${specific}_conservation.tsv.${element}; # due to memory problems, I am removing it

      echo "DONE! ${specific} ${element} ";

   done;

   # once we have used the file, we remove it
   rm ${global_accumulated};

   echo "Intersection DONE! ${specific}";
   echo "";

done;




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



# We also add PhyloCSF data as the last step

specificConservationOutput="$conservation_output_home/phyloCSF";
mkdir -p $specificConservationOutput;



element="exon"
echo "Starting phyloCSF ${element} ";

uORF_input="$batch_home/nh_m_${element}_input_$input_filename.tsv";
phylocsf_output="$specificConservationOutput/uORFs_scored_phylocsf.tsv";

python $scripts/phylocsf.py $uORF_input $phylocsf_output;



awk 'BEGIN{FS=OFS="\t"}{a[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11] = a[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11]"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19;} END{ for (i in a) print  i""a[i] }' $conservation_output_home/final_all_conservation_by_uORF.tsv $phylocsf_output > $conservation_output_home/tmp; cut -f23- --complement $conservation_output_home/tmp > $conservation_output_home/final_with_phylocsf_conservation_by_uORF.tsv; rm $conservation_output_home/tmp;



# we may modify some columns (ID) to make it easier to compare and if I could join those loci that are the same it would also be great

# 7 nucleotide length ( we must add 3 because the stop codon is not considered and we want to consider it)
# 13-21
# chr   uORFtype strand   uORF_id  gene_id  transcript_id  nuc_len  aa_len   extra cdna_seq prot_seq existingStartCodonIDs   start_phastCons100way start_phyloP100way   exon_phastCons100way exon_phyloP100way stop_phastCons100way stop_phyloP100way
# chrX	upstream_uORF	-	uORF-ENST00000612152-4;chrX:100637026-100637040:-	ENSG00000000003	ENST00000612152	12	4	.	ATGGTGAAGTTATAA	MVKL*	.	0.145	1.275	0.267	6.254	0.011	1.952

#using phyloP total
#awk 'BEGIN {FS=OFS="\t";} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13/3, $14, $15/($7+3), $16, $17/3, $18}' \
#      $filt/final_all_conservation_by_uORF.tsv > $filt/final_all_conservation_by_uORF.average.tsv;

##using phyloP average
# awk 'BEGIN {FS=OFS="\t";} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13/3, $14/3, $15/($7+3), $16/($7+3), $18/3, $19/3}' \
#      $conservation_output_home/final_all_conservation_by_uORF.tsv > $conservation_output_home/final_all_conservation_by_uORF.average.tsv;
# awk 'BEGIN {FS=OFS="\t";} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12":"$17, $13/3, $14/3, $15/($7+3), $16/($7+3), $18/3, $19/3}' \
#      $conservation_output_home/final_all_conservation_by_uORF.tsv > $conservation_output_home/final_all_conservation_by_uORF.average.tsv;
awk 'BEGIN {FS=OFS="\t";} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12":"$17, $13/3, $14/3, $15/($7+3), $16/($7+3), $18/3, $19/3, $20/($7+3), $21/($7+3), $22/($7+3)}' \
      $conservation_output_home/final_with_phylocsf_conservation_by_uORF.tsv > $conservation_output_home/final_with_phylocsf_conservation_by_uORF.average.tsv;

echo "Compressing previous files";
tar -czvf $conservation_output_home/intermediate_uORFs_conservation_files.tar.gz $conservation_output_home/all_conservation_by_* $conservation_output_home/phastCons100way/ $conservation_output_home/phyloP100way/  $conservation_output_home/unifying/ $phylocsf_output;

echo "Compressed now we remove them";
rm -rf $conservation_output_home/all_conservation_by_* $conservation_output_home/phastCons100way/ $conservation_output_home/phyloP100way/  $conservation_output_home/unifying/ $phylocsf_output;


echo "--- DONE ---";
echo "uORFs with conservation data can be found here:";
echo "$conservation_output_home/final_with_phylocsf_conservation_by_uORF.average.tsv";




