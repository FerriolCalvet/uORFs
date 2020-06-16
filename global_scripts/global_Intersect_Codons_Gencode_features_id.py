#!/usr/bin/env python

# python Intersect_StartCodons_Gencode.py /home/fcalvet/Desktop/uORFs/DATA/compact/uORFs/cleaned_compact_allENSG_25-2-2020_at_15-57.tsv /home/fcalvet/Desktop/uORFs/DATA/compact/StartCodons/cleaned_compact_reduced_start_codon.gff3 /home/fcalvet/Desktop/uORFs/DATA/merged_withStartCodons/start_codon_allENSG_25-2-2020_at_15-57.tsv




# In[0]:

import sys
import pandas as pd



# In[1]:
# print(sys.argv)
codon_type = sys.argv[1]
codon_type_split = "split_" + codon_type

# uorfs_input = '/home/fcalvet/Desktop/uORFs/DATA/compact/uORFs/cleaned_compact_allENSG_25-2-2020_at_15-57.tsv'cleaned_compact_no_hypo_only_start_codon_allENSG_23-3-2020_at_14-26.tsv
# uorfs_input = '/home/fcalvet/Desktop/uORFs/DATA/compact/uORFs/cleaned_compact_allENSG_25-2-2020_at_15-57.tsv'
uorfs_input = sys.argv[2]

# start_codons_input = '/home/fcalvet/Desktop/uORFs/DATA/compact/StartCodons/cleaned_compact_reduced_start_codon.gff3'
st_codons_input = sys.argv[3]

# merged_output = '/home/fcalvet/Desktop/uORFs/DATA/merged_withStartCodons/start_codon_allENSG_25-2-2020_at_15-57.tsv'
merged_output = sys.argv[4]



# In[2]: Read the input files and store into pandas dataframe

uorfs = pd.read_csv(uorfs_input, sep = "\t",
                      names = ['chr', 'uORFtype', 'feature_type', 'start', 'end', 'strand', 'uORF_id', 'gene_id', 'transcript_id',
                                 'nuc_len', 'aa_len', 'extra', 'cdna_seq', 'prot_seq'])
# uorfs['end'] = uorfs['end'].astype(int)

st_codons = pd.read_csv(st_codons_input, sep = "\t",
                            names = ['chr', 'source', 'start', 'end', 'strand', 
                               'existingStCodonIDs'])

# print(uorfs.shape)
# print(uorfs.size)
# uorfs.head()


# In[3]: Although they should be sorted, make sure they are

uorfs.sort_values(by = ['chr', 'start', 'end'], inplace = True)
uorfs.reset_index(drop=True, inplace = True)

st_codons.sort_values(by = ['chr', 'start', 'end'], inplace = True)
st_codons.reset_index(drop=True, inplace = True)



# In[4]: Store the identifier summarizing a start codon

st_codon_ids = {}
for index, row in st_codons.iterrows():
    st_codon_ids["{}:{}-{}:{}".format(row['chr'], row['start'],row['end'],row['strand'])] = index


# In[5]:

# print(st_codon_ids)
# uorfs.size


# In[6]: Initialize the two arrays that will be used for storing the information
#        that will be pasted to the uORFs dataframe

# not quite sure if we need this, I am not keeping it, but if needed, uncomment it
#uorfs_add_stCodon = ["." for val in range(uorfs.shape[0])]

uorfs_add_IDs = ["." for val in range(uorfs.shape[0])] 


# In[7]: Iterate over the uorfs dataframe looking for rows that have an already
#        existing start codon

st_candidates = list(st_codon_ids.keys())

for index, row in uorfs.iterrows():
    if codon_type_split not in row["extra"]:
        temp_st = "{}:{}-{}:{}".format(row['chr'],row['start'],row['end'],row['strand'])
#        uorfs_add_stCodon[index] = temp_st

        if temp_st in st_candidates:
            index_in_sts = st_codon_ids[temp_st]
            uorfs_add_IDs[index] = st_codons.iloc[index_in_sts,:]['existingStCodonIDs']
            # could be optimized using indices to refer to the field we are interested in

#            print(row['transcript_id'], temp_st)


# In[8]:

# # print(type(uorfs_add_stCodon))
# print(type(uorfs_add_IDs))


# In[9]:

host_transcripts = pd.DataFrame(uorfs_add_IDs)
all_uorfs = pd.concat([uorfs, host_transcripts], axis=1, ignore_index=True)
all_uorfs.columns = ['chr', 'uORFtype', 'feature_type', 'start', 'end', 'strand', 'uORF_id',
                     'gene_id',
                     'transcript_id', 'nuc_len', 'aa_len', 'extra', 'cdna_seq', 'prot_seq', # 'st_codon',
                     'existingStCodonIDs']


# In[10]:

all_uorfs.to_csv(merged_output, sep = "\t", header = False, index= False)






##### To merge the start codons intersected with the file containing all the uORFs information.

#awk '
#    /start_codon/{                
#        getline <"/home/fcalvet/Desktop/uORFs/DATA/merged_withStartCodons/only_start_codon_splitted_allENSG_12-3-2020_at_12-36.tsv"  #read 1 line from outer file into $0 
#    }
#    1                         #alias for `print $0`
#    ' /home/fcalvet/Desktop/uORFs/DATA/compact/uORFs/cleaned_compact_splitted_allENSG_12-3-2020_at_12-36.tsv > /home/fcalvet/Desktop/uORFs/DATA/merged_withStartCodons/all_start_codon_splitted_allENSG_12-3-2020_at_12-36.tsv

# grep -v hypothetical_start all_start_codon_splitted_allENSG_12-3-2020_at_12-36.tsv > all_start_codon_splitted_allENSG_12-3-2020_at_12-36_no_hypo_sc_separated.tsv 




