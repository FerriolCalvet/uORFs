# I have only done it for gff3 files, at some point I may do it for gtf
# gencode_file='/home/fcalvet/Desktop/REFERENCE_DATA/gencode.v33.annotation.gff3'
# gencode_file='/home/fcalvet/Desktop/REFERENCE_DATA/gencode.v33.annotation.gtf'

# It would be nice to create a variable called HOME_DIRECTORY in such a way that all
#  the structure is created inside that location, making it universal to all computers
#  now it is uncomfortable to move it

# Receive the input from the command line arguments
#      $1, $2, $3 ...
home="/home/fcalvet/Desktop/uORFs";
codon_type=$1; # either start or stop

# Create the folders needed
# mkdir -p "$home/uORFs" "$home/uORFs/SCRIPTS" "$home/uORFs/DATA";
# mkdir -p "$home/DATA/raw_data/raw_StartCodons";
# mkdir -p "$home/DATA/compact/StartCodons";

mkdir -p "$home/DATA/raw_data/raw_${codon_type}Codons";
mkdir -p "$home/DATA/compact/${codon_type}Codons";


# Define paths
#gencode_file="${home}/REFERENCE_DATA/gencode.v33.annotation.gff3";
gencode_file="${home}/REFERENCE_DATA/gencode.v34.annotation.gff3";
codons_all_info="${home}/DATA/raw_data/raw_${codon_type}Codons/${codon_type}_codon.gff3";
python_script_input="${home}/DATA/raw_data/raw_${codon_type}Codons/reduced_${codon_type}_codon.gff3";

# This command takes the gff3, subsets the start_codon lines, and saves the information we
# are interested in, in the proper format.
# It also removes gene and transcript versions
grep "${codon_type}_codon" $gencode_file > $codons_all_info;
sed -e 's/[;=]/\t/g' $codons_all_info | cut -f 1,2,4,5,7,14,16 | sed 's/\./\t/g' | cut -f 1-6,8 > $python_script_input;



python_script="${home}/SCRIPTS/global_scripts/group_Codons.py";
python_script_output="${home}/DATA/compact/${codon_type}Codons/compact_reduced_${codon_type}_codon.gff3";

# Run this python script to group the lines that are identical uORFs with different transcript/gene IDs
python $python_script $python_script_input $python_script_output;




compacting_output="${home}/DATA/compact/${codon_type}Codons/cleaned_compact_reduced_${codon_type}_codon.gff3";

# Remove all the repeated words in the same lines, keeping only the first occurrence of each
sed -r ':a; s/\b([[:alnum:]]+)\b(.*)\b\1\b/\1\2/g; ta; s/(, )+/, /g; s/, *$//' $python_script_output > $compacting_output.tmp;

# Remove the extra commas in the file by collapsing all the repeated occurrences and removing those
# commas at the end of a column (just before the tab)
# I have not checked if there can be commas at the end of the line, maybe I should
tr -s ',' < $compacting_output.tmp | sed 's/,\t/\t/g' | sed 's/,\;/\;/g' > $compacting_output;
# tr -s ',' < $compacting_output.tmp | sed 's/,\t/\t/g' > $compacting_output;

# Remove the temporary file we created to store the half cleaned version of the compact form
rm $compacting_output.tmp;


# Remove duplicated geneIDs
# not working properly (?:^|\G)(\b\w+\b),?(?=.*\1)

# Commands used to compact and clean the files after merging
# sed -r ':a; s/\b([[:alnum:]]+)\b(.*)\b\1\b/\1\2/g; ta; s/(, )+/, /g; s/, *$//' transcripts_genes_only_protein_coding.tsv > compact_transcripts_genes_only_protein_coding.tsv
# tr -s ',' < compact_transcripts_genes_only_protein_coding.tsv | sed 's/,\t/\t/g' > cleaned_compact_transcripts_genes_only_protein_coding.tsv
# https://stackoverflow.com/questions/30294915/how-to-remove-duplicate-words-from-a-string-in-a-bash-script
# https://unix.stackexchange.com/questions/187052/how-to-remove-duplicate-characters/187055



