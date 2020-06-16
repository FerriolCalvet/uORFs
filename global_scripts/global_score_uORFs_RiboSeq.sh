home="/home/fcalvet/Desktop/uORFs";


riboseq_home="${home}/REFERENCE_DATA/riboseq";
riboseq_software="/home/fcalvet/Desktop/uORFs/REFERENCE_DATA/riboseq/software";

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
riboseq_output_home="$batch_home/riboseq";
mkdir -p $riboseq_output_home;



for element in start stop exon; do
#for element in exon; do
   echo "Intersection Start ${element}";
   #   element="start";
   #   element="stop";
   #   element="exon";

   uORF_input="$batch_home/nh_m_${element}_input_$input_filename.tsv";



   # this loop allows us to process the RP coverage data and merge it with the uORFs candidates
   for specific in initiating elongating_A_site elongating_footprints; do
      echo "Starting ${specific} ${element} ";
   #   specific=initiating;
   #   specific=elongating_A_site;
   #   specific=elongating_footprints;

      if [ "$element" != "start" ] && [ "$specific" == "initiating" ]; then
         continue;
      fi

      specificRiboSeqOutput="${riboseq_output_home}/${specific}RiboSeq";
      mkdir -p $specificRiboSeqOutput;

      specificRiboSeqInput="${riboseq_home}/${specific}RiboSeq";
      global_accumulated="$specificRiboSeqInput/accumulated_${specific}.bg";



      ### START intersect

      # columns need to be adjusted to the specific cases we are handling, this is for initiating
      # intersection between the bedgraph with all the coverage of all files, and the uORF candidates
      # I flipped a and b to have a proper combination and output with the "wao" intersection
      # # bedtools intersect -wao -a $global_accumulated -b $specificRiboSeq/uORFs_input_${specific}.tsv > $specificRiboSeq/uORFs_${specific}_RPcoverage_pre.tsv
      echo "Starting bedtools intersect";
      bedtools intersect -sorted -wao -a $uORF_input -b $global_accumulated > $specificRiboSeqOutput/uORFs_${specific}_RPcoverage_pre.tsv.${element};
      echo "Finishing bedtools intersect";

      echo "Multiplying";
      awk '{ print $0"\t"$19 * $20 }' $specificRiboSeqOutput/uORFs_${specific}_RPcoverage_pre.tsv.${element} > $specificRiboSeqOutput/uORFs_${specific}_RPcoverage.tsv.${element};
      rm $specificRiboSeqOutput/uORFs_${specific}_RPcoverage_pre.tsv.${element};

      echo "Joining";
      awk 'BEGIN{FS=OFS="\t"}{a[$1"\t"$2+1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15]+=$21;}END{ for (i in a) print  i"\t"a[i] }' $specificRiboSeqOutput/uORFs_${specific}_RPcoverage.tsv.${element} | sort -k1,1 -k2,2n > $specificRiboSeqOutput/uORFs_${specific}_RPcoverage_joined.tsv.${element};
      rm $specificRiboSeqOutput/uORFs_${specific}_RPcoverage.tsv.${element};

      # process the file to recompute the total coverage of each region, we do it by multiplying the
      # coverage value of an entry times the position with overlap


      echo "DONE! ${specific} ${element} ";

   done;

   echo "Intersection DONE! ${element}";

done;



for element in start stop exon; do
   echo "Joining all ${element} in this order: elongating_A_site elongating_footprints initiating";

   if [ "$element" == "start" ]; then

      awk 'BEGIN{FS=OFS="\t"}{a[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15]=a[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15]"\t"$16;}END{ for (i in a) print  i""a[i] }' $riboseq_output_home/elongating_A_siteRiboSeq/uORFs_elongating_A_site_RPcoverage_joined.tsv.${element} $riboseq_output_home/elongating_footprintsRiboSeq/uORFs_elongating_footprints_RPcoverage_joined.tsv.${element} $riboseq_output_home/initiatingRiboSeq/uORFs_initiating_RPcoverage_joined.tsv.${element} | sort -k1,1 -k2,2n > $riboseq_output_home/all_RiboSeqcoverage_by_${element}.tsv;


   else

      # here I add a . to account for the lack of intersection with initiating coverage
      awk 'BEGIN{FS=OFS="\t"}{a[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15]=a[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15]"\t"$16;}END{ for (i in a) print  i""a[i]"\t." }' $riboseq_output_home/elongating_A_siteRiboSeq/uORFs_elongating_A_site_RPcoverage_joined.tsv.${element} $riboseq_output_home/elongating_footprintsRiboSeq/uORFs_elongating_footprints_RPcoverage_joined.tsv.${element} | sort -k1,1 -k2,2n > $riboseq_output_home/all_RiboSeqcoverage_by_${element}.tsv;

   fi

done;


# merge multiple exons into a unique entry, same for splitted start and stop codons
unifying="$riboseq_output_home/unifying";
mkdir -p $unifying;

for element in start stop exon; do
   echo "combining same feature information in ${element}";

   awk 'BEGIN{FS=OFS="\t"}{  \
         a[$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15]+=$16; \
         b[$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15]+=$17; \
         c[$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15]+=$18; \
         }END{ for (i in a) print  i"\t"a[i]"\t"b[i]"\t"c[i] }' $riboseq_output_home/all_RiboSeqcoverage_by_${element}.tsv | sort -k1,1 -k2,2n > $unifying/all_RiboSeqcoverage_unified_${element}.tsv;

done;



echo "Joining start, exon and stop in this order";
awk 'BEGIN{FS=OFS="\t"}{a[$1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12]=a[$1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12]"\t"$13"\t"$14"\t"$15"\t"$16;}END{ for (i in a) print  i""a[i] }' $unifying/all_RiboSeqcoverage_unified_start.tsv $unifying/all_RiboSeqcoverage_unified_exon.tsv $unifying/all_RiboSeqcoverage_unified_stop.tsv | sort -k5,5 -k4,4 > $riboseq_output_home/all_RiboSeqcoverage_by_uORF.tsv;



echo "Remove the column with exon matching reference";
# cut -f 16,20 --complement $riboseq_output_home/all_RiboSeqcoverage_by_uORF.tsv > $riboseq_output_home/final_all_RiboSeqcoverage_by_uORF.tsv;
cut -f 16 --complement $riboseq_output_home/all_RiboSeqcoverage_by_uORF.tsv > $riboseq_output_home/final_all_RiboSeqcoverage_by_uORF.tsv;

rm $riboseq_output_home/all_RiboSeqcoverage_by_uORF.tsv;


# 7 nucleotide length ( we must add 3 because the stop codon is not considered and we want to consider it)
# 13-21
# chr   uORFtype  strand   uORF_id  gene_id  transcript_id  nuc_len  aa_len   extra cdna_seq prot_seq existingStartCodonIDs   start_elongating_A   start_elongating_footprints   start_initiating  exon_elongating_A exon_elongating_footprints exon_initiating   stop_elongating_A stop_elongating_footprints stop_initiating

echo "Computing average";
#awk 'BEGIN {FS=OFS="\t";} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13/3, $14/3, $15/3, $16/($7+3), $17/($7+3), $18/($7+3), $19/3, $20/3, $21/3}' \
#      $filt/final_all_RiboSeqcoverage_by_uORF.tsv > $filt/final_all_RiboSeqcoverage_by_uORF.average.tsv;

awk 'BEGIN {FS=OFS="\t";} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12":"$19, $13/3, $14/3, $15/3, $16/($7+3), $17/($7+3), $18/($7+3), $20/3, $21/3, $22/3}' \
      $riboseq_output_home/final_all_RiboSeqcoverage_by_uORF.tsv > $riboseq_output_home/final_all_RiboSeqcoverage_by_uORF.average.tsv;


echo "Compressing previous files";
tar -czvf $riboseq_output_home/intermediate_uORFs_Riboseq_files.tar.gz $riboseq_output_home/all_RiboSeqcoverage_by_* $riboseq_output_home/elongating_A_siteRiboSeq/ $riboseq_output_home/elongating_footprintsRiboSeq/ $riboseq_output_home/initiatingRiboSeq/ $riboseq_output_home/unifying/;

echo "Compressed now we remove them";
rm -rf $riboseq_output_home/all_RiboSeqcoverage_by_* $riboseq_output_home/elongating_A_siteRiboSeq/ $riboseq_output_home/elongating_footprintsRiboSeq/ $riboseq_output_home/initiatingRiboSeq/ $riboseq_output_home/unifying/;

echo "--- DONE ---";
echo "uORFs with riboseq data can be found here:";
echo "$riboseq_output_home/final_all_RiboSeqcoverage_by_uORF.average.tsv";


