#!/usr/bin/env python


import sys
import pandas as pd

# sed -e 's/[;=]/\t/g' start_codon.gff3 | cut -f 1,2,4,5,7,14,16 > reduced_start_codon.gff3

# Commands used to compact and clean the files after merging
# sed -r ':a; s/\b([[:alnum:]]+)\b(.*)\b\1\b/\1\2/g; ta; s/(, )+/, /g; s/, *$//' transcripts_genes_only_protein_coding.tsv > compact_transcripts_genes_only_protein_coding.tsv
# tr -s ',' < compact_transcripts_genes_only_protein_coding.tsv | sed 's/,\t/\t/g' > cleaned_compact_transcripts_genes_only_protein_coding.tsv
# https://stackoverflow.com/questions/30294915/how-to-remove-duplicate-words-from-a-string-in-a-bash-script
# https://unix.stackexchange.com/questions/187052/how-to-remove-duplicate-characters/187055



# In[2]: Define input and output files

# input_path = "/home/fcalvet/Desktop/uORFs/DATA/raw_data/raw_uORFs/"
# input_file = "allENSG_25-2-2020_at_15-57.tsv"
# input_filename = '/home/fcalvet/Desktop/uORFs/uORF_data/withTSS/start_codon.gtf'
# input_filename = input_path + input_file

# output_path = "/home/fcalvet/Desktop/uORFs/DATA/compact/uORFs/"
# output_file = "compact_allENSG_25-2-2020_at_15-57.tsv"
# output_filename = '/home/fcalvet/Desktop/uORFs/uORF_data/withTSS/start_codons_merged.gtf'
# output_filename = output_path + output_file

input_filename = sys.argv[1]
output_filename = sys.argv[2]


# In[3]: Read the input file

start_codons = pd.read_csv(input_filename, sep = "\t",
                      names = ['chr', 'source', 'start', 'end', 'strand', 
                               'annotated_gene_id', 'annotated_transcript_id'])

# In[4]: Remove Start Codons of size != 3

for index, row in start_codons.iterrows():
    if row['end'] != row['start']+2:
        start_codons.drop(index, inplace = True)
start_codons.reset_index(drop=True, inplace = True)

# print(start_codons.shape)
# start_codons.head()


# In[5]: Build the grouped object, use all columns but the columns with IDs and source

group = start_codons.groupby(['chr','start','end','strand'])


# In[6]: Join the IDs for the start codons that are identical
#        same for the source column

start_codons['annotated_gene_id'] = group['annotated_gene_id'].transform(lambda x: ','.join(x))
start_codons['annotated_transcript_id'] = group['annotated_transcript_id'].transform(lambda x: ','.join(x))
start_codons['source'] = group['source'].transform(lambda x: ','.join(x))


# In[7]: Remove the duplicated rows after the joining process

start_codons.drop_duplicates(inplace= True)


# In[8]: Paste the two ID columns into a single one

start_codons['existingStartCodonIDs'] = start_codons['annotated_gene_id'] + ";" + start_codons['annotated_transcript_id']
start_codons.drop(labels = ['annotated_transcript_id','annotated_gene_id'], axis = 1, inplace = True)


# In[9]: Sort the rows by genomic coordinates and reset the indices

start_codons.sort_values(by = ['chr', 'start', 'end'], inplace = True)
start_codons.reset_index(drop=True, inplace = True)


# In[10]: Save to the output file

start_codons.to_csv(output_filename, sep = "\t", header = False, index= False)



# Commands used to compact and clean the files after merging
# sed -r ':a; s/\b([[:alnum:]]+)\b(.*)\b\1\b/\1\2/g; ta; s/(, )+/, /g; s/, *$//' transcripts_genes_only_protein_coding.tsv > compact_transcripts_genes_only_protein_coding.tsv
# tr -s ',' < compact_transcripts_genes_only_protein_coding.tsv | sed 's/,\t/\t/g' > cleaned_compact_transcripts_genes_only_protein_coding.tsv
# https://stackoverflow.com/questions/30294915/how-to-remove-duplicate-words-from-a-string-in-a-bash-script
# https://unix.stackexchange.com/questions/187052/how-to-remove-duplicate-characters/187055



