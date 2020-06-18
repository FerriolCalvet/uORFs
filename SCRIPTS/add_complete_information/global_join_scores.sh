home="/home/fcalvet/Desktop/uORFs";


riboseq_home="${home}/REFERENCE_DATA/riboseq";
riboseq_software="/home/fcalvet/Desktop/uORFs/REFERENCE_DATA/riboseq/software";

input_filename="$1";
filt="${home}/DATA/filtered";
batch_home="$filt/$input_filename";

all_scores="${home}/DATA/all_scores";
mkdir -p $all_scores;

all_scores_batch="${all_scores}/$input_filename.phylocsf";
mkdir -p $all_scores_batch;


filt_RiboSeq="${batch_home}/riboseq";
filt_cons="${batch_home}/conservation";


awk 'BEGIN{FS=OFS="\t"}{a[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12]=a[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12]"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21;}END{ for (i in a) print  i""a[i] }' $filt_RiboSeq/final_all_RiboSeqcoverage_by_uORF.average.tsv $filt_cons/final_with_phylocsf_conservation_by_uORF.average.tsv | sort -k5,5 -k4,4 > $all_scores_batch/all_joined_scores_by_uORF.average.tsv;

# In this case there is no need to clean the remaining columns as there are none
# cut -f 1-27  $all_scores_batch/all_joined_scores_by_uORF.average.tsv > $all_scores_batch/all_joined_scores_by_uORF_cleaned.average.tsv;

echo -e "chr\tuORFtype\tstrand\tuORF_id\tgene_id\ttranscript_id\tnuc_len\taa_len\textra\tcdna_seq\tprot_seq\texistingStart:Stop\tstart_elongating_A\tstart_elongating_footprints\tstart_initiating\texon_elongating_A\texon_elongating_footprints\texon_initiating\tstop_elongating_A\tstop_elongating_footprints\tstop_initiating\tstart_phastCons100way\tstart_phyloP100way\texon_phastCons100way\texon_phyloP100way\tstop_phastCons100way\tstop_phyloP100way\tphylocsf_1\tphylocsf_2\tphylocsf_3" > $all_scores_batch/all_scores_header_file.tsv;

cat $all_scores_batch/all_scores_header_file.tsv $all_scores_batch/all_joined_scores_by_uORF.average.tsv > $all_scores_batch/all_joined_scores_by_uORF.average.with_header.tsv;



# I still have to work on this with the ultimate goal of joining those uORFs that are identical
# Location and ID separated
awk 'BEGIN{FS=OFS="\t"}{a[$1"\t"$2"\t"$3"\t"$4"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12]=a[$1"\t"$2"\t"$3"\t"$4"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12]"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21;}END{ for (i in a) print  i""a[i] }' $filt_RiboSeq/final_all_RiboSeqcoverage_by_uORF.average.tsv $filt_cons/final_with_phylocsf_conservation_by_uORF.average.tsv | sed 's/;/\t/' | sort -k1,1 -k5,5 > $all_scores_batch/all_joined_scores_by_uORF.average.location.tsv;

# In this case there is no need to clean the remaining columns as there are none
# cut -f 1-29  $all_scores_batch/all_joined_scores_by_uORF.average.location.tsv > $all_scores_batch/all_joined_scores_by_uORF_cleaned.average.location.tsv;

echo -e "chr\tuORFtype\tstrand\tshort_uORF_id\tlocation\tuORF_id\tgene_id\ttranscript_id\tnuc_len\taa_len\textra\tcdna_seq\tprot_seq\texistingStart:Stop\tstart_elongating_A\tstart_elongating_footprints\tstart_initiating\texon_elongating_A\texon_elongating_footprints\texon_initiating\tstop_elongating_A\tstop_elongating_footprints\tstop_initiating\tstart_phastCons100way\tstart_phyloP100way\texon_phastCons100way\texon_phyloP100way\tstop_phastCons100way\tstop_phyloP100way\tphylocsf_1\tphylocsf_2\tphylocsf_3" > $all_scores_batch/all_scores_header_file.location.tsv;

cat $all_scores_batch/all_scores_header_file.location.tsv $all_scores_batch/all_joined_scores_by_uORF.average.location.tsv > $all_scores_batch/all_joined_scores_by_uORF.average.location.with_header.tsv;






echo "--- DONE ---";
echo "uORFs with all the features added can be found here:";
echo "With header";
echo "$all_scores_batch/all_joined_scores_by_uORF.average.location.with_header.tsv";
echo "";
echo "Without header";
echo "$all_scores_batch/all_joined_scores_by_uORF.average.location.tsv";

### To check that the labelling has worked
# cut -f7 $all_scores/all_joined_scores_by_uORF_cleaned.average.with_header.tsv | less;


# chr
#       chromosome

# uORFtype
#       type of uORF either upstream_uORF or overlapping_uORF

# strand
#       strand where we find the uORF either + or -

# uORF_id
#       ID of the uORF. uORF-TranscriptID-   (number of uORF in the transcript);chromosome:start-    end:      strand
#                       uORF-ENST00000612152-4;                                 chrX:      100637026-100637040:- 

# gene_id
#       Ensembl gene id

# transcript_id
#       Ensembl transcript id

# nuc_len
#       Nucleotide length, stop codon not included (real length = nuc_len + 3)

# aa_len
#       Length of the uORF peptide, stop codon not included

# extra
#       Flags adding information about the uORF. splitted_sc (splitted start codon), splitted_st (splitted stop codon)

# cdna_seq
#       Nucleotide sequence

# prot_seq
#       Protein sequence

# existingStartCodonIDs
#       Transcripts sharing start codon. EnsemblID;TranscriptID...

# start_elongating_A
#       Start codon RiboSeq elongating_A_site coverage divided by 3 (average score)

# start_elongating_footprints
#       Start codon RiboSeq elongating_footprints coverage divided by 3 (average score)

# start_initiating
#       Start codon RiboSeq initiating coverage divided by 3 (average score)

# exon_elongating_A
#       Exons RiboSeq elongating_A_site coverage divided by nuc_len+3 (average score)

# exon_elongating_footprints
#       Exons RiboSeq elongating_footprints coverage divided by nuc_len+3 (average score)

# exon_initiating
#       0, no intersection with initiating for exon or stop regions

# stop_elongating_A
#       Stop codon RiboSeq elongating_A_site coverage divided by 3 (average score)

# stop_elongating_footprints
#       Stop codon RiboSeq elongating_footprints coverage divided by 3 (average score)

# stop_initiating
#       0, no intersection with initiating for exon or stop regions

# start_phastCons100way
#       Start codon phastCons100way divided by 3 (average score)

# start_phyloP100way
#       Start codon phyloP100way divided by 3 (average score)

# exon_phastCons100way
#       Exons phastCons100way divided by nuc_len+3 (average score)

# exon_phyloP100way
#       Exons phyloP100way divided by nuc_len+3 (average score)

# stop_phastCons100way
#       Stop codon phastCons100way divided by 3 (average score)

# stop_phyloP100way
#       Stop codon phyloP100way divided by 3 (average score)

# phylocsf_0
#       phylocsf score of the uORF starting at frame 0 (+1/-1) in the start codon

# phylocsf_1
#       phylocsf score of the uORF starting at frame 1 (+2/-2) in the start codon

# phylocsf_2
#       phylocsf score of the uORF starting at frame 2 (+3/-3) in the start codon



