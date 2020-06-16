
home="/home/fcalvet/Desktop/uORFs";

scripts="$home/SCRIPTS/global_scripts";

filt="${home}/DATA/filtered";

all_scores_home="${home}/DATA/all_scores";


input_data=$(realpath $1);
# echo $input_data;

input_filename=$(basename $input_data .tsv);

batch_home="$filt/$input_filename";

mkdir -p $batch_home;



# Create initial files of the batch
create_ini_batch="0"

if [ "$create_ini_batch" == "1" ]; then
   # Remove hypothetical uORFs
   grep -v hypothetical $input_data > $batch_home/no_hypo_$input_filename.tsv;
   echo "Non hypothetical file has been created";

   # Separate start, stop and exons
   grep -w start_codon $batch_home/no_hypo_$input_filename.tsv > $batch_home/nh_start_$input_filename.tsv;
   grep -w stop_codon $batch_home/no_hypo_$input_filename.tsv > $batch_home/nh_stop_$input_filename.tsv;
   grep -w exon $batch_home/no_hypo_$input_filename.tsv > $batch_home/nh_exon_$input_filename.tsv;
   echo "Non hypothetical start, exon, stop files have been created";
fi


# Create start and stop codons reference files
create_ref="0";

if [ "$create_ref" == "1" ]; then
   bash $scripts/create_compactCodons.sh start;
   bash $scripts/create_compactCodons.sh stop;
   echo "Reference start and stop codons files have been created";
fi



# Merge start/stop codons with the reference start/stop codons from the gencode file
intersect_cod="0";

if [ "$intersect_cod" == "1" ]; then
   python_script="$scripts/global_Intersect_Codons_Gencode_features_id.py";
   start_codons_ref="$home/DATA/compact/startCodons/cleaned_compact_reduced_start_codon.gff3";
   stop_codons_ref="$home/DATA/compact/stopCodons/cleaned_compact_reduced_stop_codon.gff3";

   python $python_script \
            start $batch_home/nh_start_$input_filename.tsv \
            $start_codons_ref \
            $batch_home/nh_start_merged_$input_filename.tsv;
   echo "Start codons merged";

   python $python_script \
            stop $batch_home/nh_stop_$input_filename.tsv \
            $stop_codons_ref \
            $batch_home/nh_stop_merged_$input_filename.tsv;
   echo "Stop codons merged";

   sed "s/$/\t./" $batch_home/nh_exon_$input_filename.tsv > $batch_home/nh_exon_merged_$input_filename.tsv;
   echo "Exons merged";
fi




# we only perform this step if create is equal to 1, otherwise we assume that the input files have already been created
create_input="0";

if [ "$create_input" == "1" ]; then
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

   echo "Compressing previous files";
   tar -czvf $batch_home/initial_uORFs_files.tar.gz $batch_home/no_hypo_* $batch_home/nh_exon_* $batch_home/nh_stop_* $batch_home/nh_start_*;

   echo "Compressed now we remove them";
   rm $batch_home/no_hypo_* $batch_home/nh_exon_* $batch_home/nh_stop_* $batch_home/nh_start_*;

fi








# Add RiboSeq scores
add_RiboSeq="0";
if [ "$add_RiboSeq" == "1" ]; then
   echo "";
   echo "STARTING RIBOSEQ";
   # Add riboseq coverage features to the uORFs
   bash $scripts/global_score_uORFs_RiboSeq.sh $input_filename;
fi

# Add conservation scores
add_conservation="0";
if [ "$add_conservation" == "1" ]; then
   echo "";
   echo "STARTING CONSERVATION";
   # Add conservation features to the uORFs
   bash $scripts/global_score_uORFs_conservation.sh $input_filename;
fi

# Join RiboSeq and conservation scores
join_RiboSeq_conservation="0";
if [ "$join_RiboSeq_conservation" == "1" ]; then
   echo "";
   echo "STARTING JOINING";
   # Join conservation and riboseq coverage features
   bash $scripts/global_join_scores.sh $input_filename;
fi


# Join identical uORFs

# awk 'BEGIN{FS=OFS="\t"}{a[$1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12]=a[$1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12]"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21;}END{ for (i in a) print  i""a[i] }' $filt_RiboSeq/final_all_RiboSeqcoverage_by_uORF.average.tsv $filt_cons/final_all_conservation_by_uORF.average.tsv | sed 's/;/\t/' | sort -k1,1 -k5,5 > $all_scores_batch/all_joined_scores_by_uORF.average.location.tsv;



# Join identical uORFs


# Add additional information
join_identical_plus_features="0";
if [ "$join_identical_plus_features" == "1" ]; then
   # Run python script to join uORFs and add homology and expression data
   python $scripts/FinalJoiningStep.py \
          $all_scores_home/$input_filename/all_joined_scores_by_uORF_cleaned.average.location.with_header.tsv  \
          $all_scores_home/$input_filename/complete_uORFs_with_exp_homATG.tsv
#          $home/all_scores/$input_filename/complete_uORFs_all_near_cognate_with_exp_homATG.tsv
#          $home/all_scores/$input_filename/complete_uORFs_with_exp_hom_all_starts.tsv
#          $home/all_scores/$input_filename/complete_uORFs_all_near_cognate_with_exp_hom_all_starts.tsv
# these are the suggested filenames to keep consistency

fi




