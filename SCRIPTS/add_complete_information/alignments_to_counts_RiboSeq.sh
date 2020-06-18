
home="/home/fcalvet/Desktop/uORFs";

riboseq_home="${home}/REFERENCE_DATA/riboseq";
riboseq_software="/home/fcalvet/Desktop/uORFs/REFERENCE_DATA/riboseq/software";



# $riboseq_software/bigWigToBedGraph /home/fcalvet/Desktop/uORFs/REFERENCE_DATA/conservation/hg38.phastCons100way.bw /home/fcalvet/Desktop/uORFs/REFERENCE_DATA/conservation/hg38.phastCons100way.bg

# $riboseq_software/bigWigToBedGraph /home/fcalvet/Desktop/uORFs/REFERENCE_DATA/conservation/hg38.phyloP100way.bw /home/fcalvet/Desktop/uORFs/REFERENCE_DATA/conservation/hg38.phyloP100way.bg



# this loop is to change the file format of the files containing coverage data
# from bigWig to BedGraph

for specific in initiating elongating_A_site elongating_footprints; do
   specificRiboSeq="${riboseq_home}/${specific}RiboSeq";
   # change the file format from bigWig to BedGraph
   for var in $(ls $specificRiboSeq/bw | cut -d '.' -f1,2); do
      $riboseq_software/bigWigToBedGraph $specificRiboSeq/bw/$var.bw $specificRiboSeq/bg/$var.bg;
   done;
done;



# this loop allows us to process the RP coverage data and merge it with the uORFs candidates
for specific in initiating elongating_A_site elongating_footprints; do

#   specific=initiating;
#   specific=elongating_A_site;
#   specific=elongating_footprints;

   specificRiboSeq="${riboseq_home}/${specific}RiboSeq"; 


   ### START join all the data from the different files into a single file

   # Create the string with all the files for the command
   store="";
   bg_path="$specificRiboSeq/bg";
   output_file="$specificRiboSeq/bg_joined_${specific}.bg";

   for new_var in $(ls $bg_path); do
      n_store="$store $bg_path/$new_var";
      store=$n_store;
   done;

   # Run the command, it will merge all the bedgraph files with coverage data into a single one
   # with the coordinates in the first 3 columns, and then as many columns as input files
   bedtools unionbedg -i $store > $output_file;


   # Add all the lines providing information about the counts and create a new file containing only the positions and the sum
   # It will still be in BedGraph format
   global_accumulated="$specificRiboSeq/accumulated_${specific}.bg";

            # Source
            # http://www.unixcl.com/2009/12/sum-numbers-in-each-row-awk.html
   awk 'BEGIN {FS=OFS="\t"} {sum=0; for(i=4;i<=NF;i++) {sum+=$i;} print $1, $2, $3, sum}' $output_file > $global_accumulated;

   ### END join all the data from the different files into a single file

   rm $output_file;


done;

