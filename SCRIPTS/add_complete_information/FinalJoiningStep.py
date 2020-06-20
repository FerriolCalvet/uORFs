#!/usr/bin/env python
# coding: utf-8

import sys, os
import pandas as pd


# In[1]: Define input and output files

# input_filename = '/home/fcalvet/Desktop/uORFs/DATA/all_scores/allENSG_20-4-20_at_14-22/all_joined_scores_by_uORF_cleaned.average.location.with_header.tsv'
input_filename = os.path.abspath(sys.argv[1])
short_name = input_filename.split(sep = '/')[-1].split(sep = '.')[0]
print(short_name)


# output_filename = '/home/fcalvet/Desktop/uORFs/DATA/all_scores/allENSG_20-4-20_at_14-22/complete_uORFs_with_expression_mouse_homology.tsv'
output_filename = sys.argv[2]

# Reference files
# home="/home/fcalvet/Desktop/uORFs"
home="../.."

ref_home = home+"/REFERENCE_DATA"

# expression_file="/home/fcalvet/Desktop/yale_pipeline/GTEX_tissue_meanexpression_entropy_no_version.tsv"
expression_file=ref_home+"/expression/GTEX_tissue_meanexpression_entropy_no_version.tsv"

homologs_list=ref_home+'/orthologues/human_to_mouse_ortholog_one2one_high_conf.tsv'

genes_with_uorfs_mouse=ref_home+'/orthologues/genes_with_uORFs_in_mouse'

# For all near cognate start codons
#genes_with_uorfs_mouse=ref_home+'/orthologues/genes_with_uORFs_in_mouse_all_near_cognate'





# In[2]: Read the input file
uorfs_initial = pd.read_csv(input_filename, sep = "\t", index_col = False)

uorfs_initial.columns = ['chr', 'uORFtype', 'strand', 'short_uORF_id', 'location', 'uORF_id',
       'gene_id', 'transcript_id', 'nuc_len', 'aa_len', 'extra', 'cdna_seq',
       'prot_seq', 'existingStart:Stop', 'start_elongating_A',
       'start_elongating_footprints', 'start_initiating', 'exon_elongating_A',
       'exon_elongating_footprints', 'exon_initiating', 'stop_elongating_A',
       'stop_elongating_footprints', 'stop_initiating',
       'start_phastCons100way', 'start_phyloP100way', 'exon_phastCons100way',
       'exon_phyloP100way', 'stop_phastCons100way', 'stop_phyloP100way',
       'phylocsf_1', 'phylocsf_2', 'phylocsf_3']

expression = pd.read_csv(expression_file, sep = "\t", index_col = False)
# Remove version from this file
# sed 's/\(ENSG[0-9]*\)\.[0-9]*/\1/g' new_u | sed 's/\(ENST[0-9]*\)\.[0-9]*/\1/g' > c_new_u

expression.columns = ['transcript_id', 'Thyroid', 'Testis', 'Cervix Uteri', 'Adipose Tissue',
       'Breast', 'Vagina', 'Nerve', 'Pituitary', 'Stomach', 'Fallopian Tube',
       'Bone Marrow', 'Bladder', 'Blood', 'Colon', 'Prostate', 'Pancreas',
       'Blood Vessel', 'Liver', 'Spleen', 'Small Intestine', 'Uterus', 'Ovary',
       'Muscle', 'Heart', 'Adrenal Gland', 'Brain', 'Salivary Gland', 'Lung',
       'Skin', 'Esophagus', 'Kidney', 'Tissue Entropy']


homologs = pd.read_csv(homologs_list, sep = "\t", index_col = False, header = None)

homologs.columns = ['gene_id', 'mouse_gene_id', 'confidence', 'homology_type']


uorfs_mouse = pd.read_csv(genes_with_uorfs_mouse, sep = "\t", index_col = False, header = None)
uorfs_mouse.columns = ['mouse_gene_id']
uorfs_mouse['num_homologs_with_uORFs'] = 1


print('Reading done!')


# In[3]: Merge with expression and homology data

expression_added = uorfs_initial.merge(expression, on = "transcript_id" , how = "left")

with_homology_and_expression = expression_added.merge(homologs, on = "gene_id" , how = "left")

expression_added = with_homology_and_expression.merge(uorfs_mouse, on = "mouse_gene_id" , how = "left")

print('Merging done!')

del uorfs_initial
del with_homology_and_expression

# # In[4]: Save the dataframe to a file
#expression_output_filename = '/home/fcalvet/Desktop/uORFs/DATA/all_scores/all_near_cognate_ENSG_7-4-2020_at_11-54/uORFs_with_expression'
#expression_output_filename = '/home/fcalvet/Desktop/uORFs/DATA/all_scores/allENSG_23-3-2020_at_18-29/uORFs_with_expression'
#expression_added.to_csv(expression_output_filename, sep = "\t", index= False)


# In[5]: Correct the columns that have NaN values due to the lack of homologs
expression_added['num_homologs_with_uORFs'] = expression_added['num_homologs_with_uORFs'].where(pd.notnull(expression_added['num_homologs_with_uORFs']), 0)
expression_added['num_homologs_with_uORFs'] = expression_added['num_homologs_with_uORFs'].astype(int)

expression_added['confidence'] = expression_added['confidence'].where(pd.notnull(expression_added['confidence']), 0)
expression_added['confidence'] = expression_added['confidence'].astype(int)

expression_added[['mouse_gene_id','homology_type']] = expression_added[['mouse_gene_id','homology_type']].where(pd.notnull(expression_added[['mouse_gene_id', 'homology_type']]), None)


print('Homology columns corrected!')

# In[6]: This first grouping is to join the IDs and string fields

expression_grouped = expression_added.groupby(['chr', 'uORFtype', 'strand', 'location',
                                               'nuc_len', 'aa_len', 'extra', 'cdna_seq',
                                               'prot_seq', 'existingStart:Stop', 'start_elongating_A',
                                               'start_elongating_footprints', 'start_initiating', 'exon_elongating_A',
                                               'exon_elongating_footprints', 'exon_initiating', 'stop_elongating_A',
                                               'stop_elongating_footprints', 'stop_initiating',
                                               'start_phastCons100way', 'start_phyloP100way', 'exon_phastCons100way',
                                               'exon_phyloP100way', 'stop_phastCons100way', 'stop_phyloP100way',
                                               'phylocsf_1', 'phylocsf_2', 'phylocsf_3'])


# In[7]: Join the text fields for the identical uORFs

expression_added['short_uORF_id'] = expression_grouped['short_uORF_id'].transform(lambda x: ';'.join(x))

expression_added['gene_id'] = expression_grouped['gene_id'].transform(lambda x: ';'.join(x))

expression_added['transcript_id'] = expression_grouped['transcript_id'].transform(lambda x: ';'.join(x))

expression_added['mouse_gene_id'] = expression_grouped['mouse_gene_id'].transform(lambda x: ';'.join(filter(None,x)))

expression_added['homology_type'] = expression_grouped['homology_type'].transform(lambda x: ';'.join(filter(None,x)))

print('Text fields joined!')

# In[8]:

#expression_added.to_csv('/home/fcalvet/Desktop/uORFs/DATA/all_scores/allENSG_23-3-2020_at_18-29/uORFs_with_expression_ids_merged', sep = "\t", index= False)


# In[9]:

#expression_added = pd.read_csv('/home/fcalvet/Desktop/uORFs/DATA/all_scores/allENSG_23-3-2020_at_18-29/uORFs_with_expression_ids_merged', sep = "\t", index_col = False)


# In[10]: Join numerical fields
new_grouped = expression_added.groupby(['chr', 'uORFtype', 'strand', 'short_uORF_id', 'location',
       'gene_id', 'transcript_id', 'nuc_len', 'aa_len', 'extra', 'cdna_seq',
       'prot_seq', 'existingStart:Stop', 'start_elongating_A',
       'start_elongating_footprints', 'start_initiating', 'exon_elongating_A',
       'exon_elongating_footprints', 'exon_initiating', 'stop_elongating_A',
       'stop_elongating_footprints', 'stop_initiating',
       'start_phastCons100way', 'start_phyloP100way', 'exon_phastCons100way',
       'exon_phyloP100way', 'stop_phastCons100way', 'stop_phyloP100way',
       'phylocsf_1', 'phylocsf_2', 'phylocsf_3',
       'mouse_gene_id', 'homology_type'])


# In[11]: Compute the mean of the rows that need to be joined

new_expression_added = new_grouped.mean()


# In[12]: Reset the indeces to reintroduce the fields of the group by to the data frame

new_expression_added.reset_index(inplace = True, drop = False)

print('Numerical fields corrected!')

# In[13]: Correct text fields to avoid repetitions

splitted_genes = new_expression_added['gene_id'].str.split(';')
corrected_genes = splitted_genes.apply(lambda x: ','.join(sorted(list(set(filter(None, x))))))
new_expression_added['gene_id'] = corrected_genes


splitted_transcripts = new_expression_added['transcript_id'].str.split(';')
corrected_transcripts = splitted_transcripts.apply(lambda x: ','.join(sorted(list(set(filter(None, x))))))
new_expression_added['transcript_id'] = corrected_transcripts


splitted_mouse_genes = new_expression_added['mouse_gene_id'].str.split(';')
corrected_mouse_genes = splitted_mouse_genes.apply(lambda x: ','.join(sorted(list(set(filter(None, x))))))
new_expression_added['mouse_gene_id'] = corrected_mouse_genes


splitted_homology = new_expression_added['homology_type'].str.split(';')
corrected_homology = splitted_homology.apply(lambda x: ','.join(sorted(list(set(filter(None, x))))))
new_expression_added['homology_type'] = corrected_homology


# In[14]: Correct the short IDs field and extend the change to the uORFs ID field
def dividing(x):
    '''
    this function is specific for joining the uORFs IDs
    '''
    all_ids_together = ','.join(sorted(list(filter(None, x))))
    key = all_ids_together[:9]
    final = key + all_ids_together.replace(key, "")    
    return final

short_IDs = new_expression_added['short_uORF_id'].str.split(';')
corrected_short_IDs_plus_start = short_IDs.apply(dividing)
new_expression_added['short_uORF_id'] = corrected_short_IDs_plus_start
new_expression_added['uORF_id'] = new_expression_added['short_uORF_id'] + ";" + new_expression_added['location']


# In[15]:

def start_codon(x):
    # this function is to extract the start codon from the short uORF id field
    startcod = x[5:8]
    return startcod

new_expression_added['start_codon'] = new_expression_added['short_uORF_id'].apply(start_codon)

print('Text fields corrected!')


# In[16]: Use the location field to define the start and end of each uORF

new_expression_added[['chr2','position','strand2']] = new_expression_added['location'].str.split(':',expand=True)
new_expression_added[['start', 'end']] = new_expression_added['position'].str.split('-',expand=True)
new_expression_added['start'] = new_expression_added['start'].astype(int)
new_expression_added['end'] = new_expression_added['end'].astype(int)
new_expression_added.drop(['chr2', 'strand2', 'position'], axis=1, inplace = True)

print('Start and end positions added!')


# In[17]: Remove duplicates, if any, and sort the uORFs according to its genomic coordinates

new_expression_added.drop_duplicates(inplace= True)
new_expression_added.sort_values(by = ['chr','start', 'end', 'uORF_id'], inplace = True)
new_expression_added.reset_index(drop=True, inplace = True)

print('Duplicates removed and sorted!')


# In[18]: Rearrange the columns in the desired order

last_df = new_expression_added[['chr', 'start', 'end', 'strand', 'uORFtype', 'start_codon', 
                                'short_uORF_id', 'location', 'uORF_id',
                                'gene_id', 'transcript_id', 'nuc_len',
                                'aa_len', 'extra', 'cdna_seq', 'prot_seq',
                                'existingStart:Stop', 'start_elongating_A',
                                'start_elongating_footprints', 'start_initiating', 'exon_elongating_A',
                                'exon_elongating_footprints', 'exon_initiating', 'stop_elongating_A',
                                'stop_elongating_footprints', 'stop_initiating',
                                'start_phastCons100way', 'start_phyloP100way', 'exon_phastCons100way',
                                'exon_phyloP100way', 'stop_phastCons100way', 'stop_phyloP100way',
                                'phylocsf_1', 'phylocsf_2', 'phylocsf_3',
                                'mouse_gene_id', 'homology_type', 'confidence', 'num_homologs_with_uORFs',
                                'Thyroid', 'Testis', 'Cervix Uteri', 'Adipose Tissue', 'Breast',
                                'Vagina', 'Nerve', 'Pituitary', 'Stomach', 'Fallopian Tube',
                                'Bone Marrow', 'Bladder', 'Blood', 'Colon', 'Prostate', 'Pancreas',
                                'Blood Vessel', 'Liver', 'Spleen', 'Small Intestine', 'Uterus', 'Ovary',
                                'Muscle', 'Heart', 'Adrenal Gland', 'Brain', 'Salivary Gland', 'Lung',
                                'Skin', 'Esophagus', 'Kidney', 'Tissue Entropy']]
print('Rearranged and saving it!')

# In[19]: Save the dataframe to a file

#last_df.to_csv('/home/fcalvet/Desktop/uORFs/DATA/all_scores/allENSG_23-3-2020_at_18-29/complete_file_start_codon_ENS', sep = "\t", header = True, index= False)
last_df.to_csv(output_filename, sep = "\t", header = True, index= False)

print('DONE!')
print('File available at: ', output_filename)

