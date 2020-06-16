#!/usr/bin/env python
# coding: utf-8

# In[1]: Import modules and receive filename from the command args

import sys
import pandas as pd
import numpy as np
import pyBigWig as pw


filename = sys.argv[1]
output_filename = sys.argv[2]



# In[2]: Read input file, we only read the exons because they contain the whole sequence

uorfs = pd.read_csv(filename, sep = "\t", index_col = False,
                            names = ['chr', 'start', 'end', 'uORFtype', 'feature_type', 'strand',
                                     'uORF_id', 'gene_id', 'transcript_id', 'nuc_len', 'aa_len',
                                     'extra', 'cdna_seq', 'prot_seq', 'existingStart:Stop'] )



# In[4]: We compute the change in the reading frame caused by the CDS traversing the whole exon.
#           This is computed by getting the residue of the exon length divided by 3

change_between_in_out_frames = (uorfs.end.astype(int) - uorfs.start.astype(int)) % 3
uorfs["change_frames"] = change_between_in_out_frames



### Add frames to the positive strand uORFs

# we create a list of 0s assuming that all exons start at frame 0, we will update it accordingly
uorfs_exon_frame = [0 for val in range(uorfs.shape[0])]

# as we first compute the frames of the exons in the positive strand, we sort by the coordinates
#     in ascending order so that we have the exons of the same uORF_id together and properly sorted
uorfs.sort_values(by = ['uORF_id', 'chr', 'start', 'end'], inplace = True, ascending= True)
uorfs.reset_index(drop=True, inplace = True)


# we compute the frame of each exon, if they are the first exon with that ID we leave them at 0
# otherwise we update it with the change of frame caused by the previous exon
prev_id = ""
prev_frame = 0
for index, row in uorfs.iterrows():
    if row["strand"] == "+":
        if row["uORF_id"] != prev_id:
            prev_id = row["uORF_id"]
#            prev_frame = (0 + row["change_frames"]) % 3
            first_exon_frame = row["start"] % 3 # we compute the frame at which the uORF starts, and we update the value
            uorfs_exon_frame[index] = first_exon_frame
            prev_frame = (first_exon_frame + row["change_frames"]) % 3
        else:
            uorfs_exon_frame[index] = prev_frame
            prev_frame = (prev_frame + row["change_frames"]) % 3

# we assign the updated list of frames to the column in the dataframe
uorfs["exon_frame"] = uorfs_exon_frame



### Add frames to the negative strand uORFs

# we sort the uORFs in decreasing order so that we can work properly with the uORFs in the negative strand
uorfs.sort_values(by = ['uORF_id', 'chr', 'start', 'end'], inplace = True, ascending= False)
uorfs.reset_index(drop=True, inplace = True)

# as we have resorted the dataframe, we also have to update the list of exon frames so that we don't lose
#     any information, and we can update the values for the negative strand without changing those of the positive
uorfs_exon_frame = uorfs["exon_frame"]

# we perfom the same procedure as we did for the positive strand uORFs but with the negative
prev_id = ""
prev_frame = 0
for index, row in uorfs.iterrows():
    if row["strand"] == "-":
        if row["uORF_id"] != prev_id:
            prev_id = row["uORF_id"]
#            prev_frame = (0 + row["change_frames"]) % 3
            first_exon_frame = row["end"] % 3 # we compute the frame at which the uORF starts, and we update the value
            uorfs_exon_frame[index] = first_exon_frame
            prev_frame = (first_exon_frame + row["change_frames"]) % 3
        else:
            uorfs_exon_frame[index] = prev_frame
            prev_frame = (prev_frame + row["change_frames"]) % 3

# we assign the updated list of uORFs exon frames to the dataframe
uorfs["exon_frame"] = uorfs_exon_frame



####
# Compute PhyloCSF scores for each location, and store them
####

# as we will use location for the next steps, we define a column with it
uorfs["location"] = uorfs['chr'] + ":" + uorfs['start'].astype(str) + "-" + uorfs['end'].astype(str) + ":" + uorfs["strand"]


# def predict_coding(vec):
#    for v in vec:
#        if not v: continue
#        if v > 0: return "yes"
#    return "no"


regs = { "+" : [], "-" : [] }
chrom={}
starts={}
ends={}


# Although we already have the file in a pandas, we read again the file creating dictionaries that
#     will be useful for getting the phylocsf scores

for line in open(filename):
    if not line.startswith("chr"):
        continue
    fields = line.strip().split()
    (chr, start, end, strand) = (fields[0], fields[1], fields[2], fields[5])
    name = chr+":"+start+"-"+end
    chrom[name]=chr
    starts[name]=int(start)
    ends[name]=int(end)
    regs[strand].append(name)



# In[19]: We fetch the phyloCSF scores from the files either local or in the network

#rpathbase = "https://data.broadinstitute.org/compbio1/PhyloCSFtracks/hg19/latest/PhyloCSF"
#rpathbase = "https://data.broadinstitute.org/compbio1/PhyloCSFtracks/hg38/20170117/PhyloCSF"
# this is the path to the folder with PhyloCSF files, this path also includes the common part of the filename
rpathbase = "/home/fcalvet/Desktop/uORFs/REFERENCE_DATA/conservation/PhyloCSF/PhyloCSF"

scores = {}
positive_regs = list(set(regs["+"]))

scores = {}
negative_regs = list(set(regs["-"]))

for rf in ["+1","+2","+3"]:
    sys.stderr.write("Reading frame " + rf + "\n")
    rpath = rpathbase + rf + ".bw"
    bw = pw.open(rpath)
    frame_score = {}
    count = 0
    for r in positive_regs:
        count += 1
#        if(count % 500 == 0): sys.stderr.write('\tProcessed ' + str(count) + "\n")
#        sys.stderr.flush()
        score = sum(bw.values(chrom[r], starts[r], ends[r]))
#        print(chrom[r], starts[r], ends[r], score)
#        score = bw.stats(chrom[r], starts[r], ends[r])[0]
#        score = bw.stats(chrom[r], starts[r], ends[r], exact = True)[0]
        if np.isnan(score):
            score = -np.inf
        frame_score[r+':+'] = score
    scores[rf] = frame_score


for rf in ["-1","-2","-3"]:
    sys.stderr.write("Reading frame " + rf + "\n")
    rpath = rpathbase + rf + ".bw"
    bw = pw.open(rpath)
    frame_score = {}
    count = 0
    for r in negative_regs:
        count += 1
#        if(count % 500 == 0): sys.stderr.write('\tProcessed ' + str(count) + "\n")
#        sys.stderr.flush()
        score = sum(bw.values(chrom[r], starts[r], ends[r]))
#        print(chrom[r], starts[r], ends[r], score)
#        score = bw.stats(chrom[r], starts[r], ends[r])[0]
#        score = bw.stats(chrom[r], starts[r], ends[r], exact = True)[0]
        if np.isnan(score):
            score = -np.inf
        frame_score[r+':-'] = score
    scores[rf] = frame_score


# In[21]:

scores_frame = pd.DataFrame.from_dict(scores)
scores_frame.reset_index(inplace = True, drop = False)
scores_frame.columns = ['location', '+1', '+2', '+3', '-1', '-2', '-3']
scores_frame["location"] = scores_frame["location"].astype(str)
uorfs["location"] = uorfs["location"].astype(str)


# In[22]:

uorfs.merge(scores_frame, on = "location")




# ### Add phyloCSF scores to the positive strand uORFs starting each uORF at 0, 1 or 2


scores_in_frame = {}


# In[25]:

uorfs.sort_values(by = ['uORF_id', 'chr', 'start', 'end'], inplace = True, ascending= True)
uorfs.reset_index(drop=True, inplace = True)


# In[27]:


prev_id = ""
prev_frame = 0
score1,score2,score3 = 0,0,0
for index, row in uorfs.iterrows():
    if row["strand"] == "+":
        if row["uORF_id"] != prev_id:
            if prev_id != "":
                scores_in_frame[prev_id] = [score1,score2,score3]
            
            prev_id = row["uORF_id"]
            
            score1,score2,score3 = 0,0,0
            
            if row["exon_frame"] == 0:
                score1 += scores["+1"][row["location"]]
                score2 += scores["+2"][row["location"]]
                score3 += scores["+3"][row["location"]]
                
            elif row["exon_frame"] == 1:
                score1 += scores["+2"][row["location"]]
                score2 += scores["+3"][row["location"]]
                score3 += scores["+1"][row["location"]]
                
            elif row["exon_frame"] == 2:
                score1 += scores["+3"][row["location"]]
                score2 += scores["+1"][row["location"]]
                score3 += scores["+2"][row["location"]]
                
                
                
        else:

            if row["exon_frame"] == 0:
                score1 += scores["+1"][row["location"]]
                score2 += scores["+2"][row["location"]]
                score3 += scores["+3"][row["location"]]
                
            elif row["exon_frame"] == 1:
                score1 += scores["+2"][row["location"]]
                score2 += scores["+3"][row["location"]]
                score3 += scores["+1"][row["location"]]
                
            elif row["exon_frame"] == 2:
                score1 += scores["+3"][row["location"]]
                score2 += scores["+1"][row["location"]]
                score3 += scores["+2"][row["location"]]

scores_in_frame[prev_id] = [score1,score2,score3]



# NEGATIVE

# In[28]:


uorfs.sort_values(by = ['uORF_id', 'chr', 'start', 'end'], inplace = True, ascending= False)
uorfs.reset_index(drop=True, inplace = True)


# In[30]:

prev_id = ""
prev_frame = 0
score1,score2,score3 = 0,0,0
for index, row in uorfs.iterrows():
    if row["strand"] == "-":
        if row["uORF_id"] != prev_id:
            if prev_id != "":
                scores_in_frame[prev_id] = [score1,score2,score3]
            
            prev_id = row["uORF_id"]
            
            score1,score2,score3 = 0,0,0
            
            if row["exon_frame"] == 0:
                score1 += scores["-1"][row["location"]]
                score2 += scores["-2"][row["location"]]
                score3 += scores["-3"][row["location"]]
                
            elif row["exon_frame"] == 1:
                score1 += scores["-2"][row["location"]]
                score2 += scores["-3"][row["location"]]
                score3 += scores["-1"][row["location"]]
                
            elif row["exon_frame"] == 2:
                score1 += scores["-3"][row["location"]]
                score2 += scores["-1"][row["location"]]
                score3 += scores["-2"][row["location"]]
                
                
                
        else:

            if row["exon_frame"] == 0:
                score1 += scores["-1"][row["location"]]
                score2 += scores["-2"][row["location"]]
                score3 += scores["-3"][row["location"]]
                
            elif row["exon_frame"] == 1:
                score1 += scores["-2"][row["location"]]
                score2 += scores["-3"][row["location"]]
                score3 += scores["-1"][row["location"]]
                
            elif row["exon_frame"] == 2:
                score1 += scores["-3"][row["location"]]
                score2 += scores["-1"][row["location"]]
                score3 += scores["-2"][row["location"]]

scores_in_frame[prev_id] = [score1,score2,score3]                
                
#         if index > 800:
#             print(index)


# In[31]:


ids_scored = pd.DataFrame.from_dict(scores_in_frame, orient = "index")


# In[32]:

ids_scored.reset_index(inplace= True, drop = False)
#ids_scored.columns = ['uORF_id', 0, 1, 2]
ids_scored.columns = ['uORF_id', 1, 2, 3]


# In[33]:

id_grouped_df = uorfs.groupby("uORF_id").first()


# In[35]:


id_grouped_df.reset_index(inplace= True, drop = False)
id_grouped_df.drop(labels = ['start','end','feature_type','existingStart:Stop', 'change_frames', 'exon_frame', 'location'],
                  axis=1, inplace = True)


# In[36]:

scores_added = id_grouped_df.merge(ids_scored, on = "uORF_id" , how = "left")


# In[41]:

final_scores_reorganized = scores_added[['chr','uORFtype','strand','uORF_id','gene_id','transcript_id',
                                         'nuc_len','aa_len','extra','cdna_seq','prot_seq',0,1,2]]


# In[42]:

final_scores_reorganized.to_csv(output_filename, sep = "\t", header = False, index= False)

