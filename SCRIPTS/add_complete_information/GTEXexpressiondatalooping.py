#Requires GTEx expression data as input.
# Usage:
# python GTEXexpressiondatalooping.py > ../../REFERENCE_DATA/expression/GTEX_tissue_meanexpression_entropy_no_version.tsv



from math import log

import os


# home = "/home/fcalvet/Desktop"
home = "../.."

expression_home = home + "/REFERENCE_DATA/expression"

os.system("mkdir -p " + expression_home)

os.system("wget -P " + expression_home + " https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")

os.system("wget -P " + expression_home + " https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct.gz")
os.system("gunzip " + expression_home + "/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct.gz")

os.system("sed -e '1,2d' " + expression_home + "/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct > " + expression_home + "/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count_correct.gct")




labels_file = expression_home + "/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"

with open(labels_file) as tissue_legend_tsv:
	tissue_legend_dict = {}
	tissue_legend_tsv.readline() # skip the header
	current_line = tissue_legend_tsv.readline()
	while current_line:
		fields = current_line.split('\t')
		samplelabel = fields[0]
		tissue = fields[5]
		if tissue in tissue_legend_dict:
			tissue_legend_dict[tissue].append(samplelabel)
		else:
			tissue_legend_dict[tissue] = [samplelabel]
		current_line = tissue_legend_tsv.readline()

all_tissue_legend_keys = tissue_legend_dict.keys()



# I have removed the first two lines of the file downloaded from 

expression_data_file = expression_home + "/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count_correct.gct"


#Getting the entry locations associated with each key
with open(expression_data_file) as gtex_expression_data:
	current_line = gtex_expression_data.readline()
	# print current_line 
	column_labels = current_line.split()
#	print len(column_labels)
	tissue_locations_dict = {}
	for key in all_tissue_legend_keys:
#		print key
		all_patient_samples = tissue_legend_dict[key]
#		print len(all_patient_samples)
		patient_sample_locations_tissue = []
		for patient_sample in all_patient_samples:
			if patient_sample in column_labels:
				patient_sample_location = column_labels.index(patient_sample)
				patient_sample_locations_tissue.append(patient_sample_location)

		#if len(patient_sample_locations_tissue) > 0: # remove tissues without any sample # it seems that it is not working
		tissue_locations_dict[key] = patient_sample_locations_tissue





#Averaging the data associated with each key:
with open(expression_data_file) as gtex_expression_data:
	#this dict structure will contain the data
	print '\t'.join(["TranscriptID", '\t'.join(map(str,all_tissue_legend_keys)), "Tissue Entropy"])
	#skip the first line, containing the header
	gtex_expression_data.readline()
	#and the second line
	current_line = gtex_expression_data.readline()
	while current_line:
		tissue_means_list = []
		tissue_contributions_entropy = []
		if current_line.startswith('transcript_id'): # check that the header has been removed
			current_line = gtex_expression_data.readline()
		if not current_line:
			break
		for key in tissue_locations_dict:
			fields = current_line.split()
			transcript_name = fields[0]
#			print fields
#			print key, tissue_locations_dict[key]
			expression_float = [float(fields[x]) for x in tissue_locations_dict[key]]
			try:
				tissue_mean = (sum(expression_float)/len(expression_float))
			except ZeroDivisionError:
				tissue_mean = 0
			tissue_means_list.append(tissue_mean)			
		try:
			relative_expression = [float(x/sum(tissue_means_list)) for x in tissue_means_list]
			for x in relative_expression:
				if x == 0:
					tissue_contribution_entropy = 0
				else:
					tissue_contribution_entropy = -x*log(x, 2)
				tissue_contributions_entropy.append(tissue_contribution_entropy)
			tissue_entropy = sum(tissue_contributions_entropy)
		except ZeroDivisionError:
			relative_expression = float('Inf')
			tissue_entropy = 0
		print '\t'.join([transcript_name, '\t'.join(map(str,tissue_means_list)), str(tissue_entropy)])
		current_line = gtex_expression_data.readline()







# To check the error in the last line
# tail -2 Desktop/yale_pipeline/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count_correct.gct | cut -f 16000- | less -S
