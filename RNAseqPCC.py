import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools

#############################################ARGUMENTS###############################################
'''
#Testing files
file_1 = "/Users/michaeltellis/Downloads/CTSP-DLBCL1_2/Xena_matrices/CTSP-DLBCL1.star_tpm.tsv"
file_2 = "/Users/michaeltellis/Downloads/CTSP-DLBCL1_2/Xena_matrices/CTSP-DLBCL1.star_fpkm-uq.tsv"
file_3 = "/Users/michaeltellis/Downloads/CTSP-DLBCL1_2/Xena_matrices/CTSP-DLBCL1.star_fpkm.tsv"
file_4 = "/Users/michaeltellis/Downloads/CTSP-DLBCL1_2/Xena_matrices/CTSP-DLBCL1.star_counts.tsv"
'''
file_1 = sys.argv[1]
file_2 = sys.argv[2]
file_3 = sys.argv[3]
file_4 = sys.argv[4]

if len(sys.argv) != 5:
	print("incorrect arguments!")
	sys.exit(0)


############################################GLOBAL###################################################
file_name = [file_1, file_2, file_3, file_4] 
files = [None] * len(file_name)
special_case = None
combinations= []

############################################constants################################################
keywords = ['counts.tsv', 'tpm.tsv', 'fpkm.tsv', 'fpkm-uq.tsv']


def specialCase(files):
	counts = None
	for i in range(len(files)):
		if files[i].endswith(keywords[0]):
			counts = i
			break
	return counts



def getSamples(xena_file):
    file = open(xena_file, "r+")
    l = file.readline()
    file.close()
    l = l.split("\t")
    l.pop(0)
    sampleList = []
    for i in l:
        sampleList.append(i.strip())
 
    return(sampleList)

def formatDf(file):
	file_df = pd.read_csv(file, sep= '\t', index_col = 0)

	return file_df


def PCCcol(df_1, df_2, sample_list, file_1, file_2):
	corr_1 = np.exp2(df_1[sample_list[0]]) - 1
	corr_2 = np.exp2(df_2[sample_list[0]]) - 1
	corr_value = corr_1.corr(corr_2)
	corr_value_non_convert = df_1[sample_list[0]].corr(df_2[sample_list[0]], method = "spearman")
	file_1_name = None
	file_2_name = None

	for key in keywords:
		if file_1.endswith(key):
			file_1_name = key
		if file_2.endswith(key):
			file_2_name = key


	print(file_1_name)
	print(file_2_name)
	print(corr_value)
	print(corr_value_non_convert)


	fig , plot = plt.subplots()
	plot.scatter(corr_1, corr_2, s = 5)
	plt.xlabel(file_1_name, fontsize = 10)
	plt.ylabel(file_2_name, fontsize = 10)
	plt.show()

def PCCrow(df_1, df_2, file_1, file_2):
	corr_1 = np.exp2(df_1.iloc[0]) - 1
	corr_2 = np.exp2(df_2.iloc[0]) - 1
	corr_value = corr_1.corr(corr_2)
	file_1_name = None
	file_2_name = None
	for key in keywords:
		if file_1.endswith(key):
			file_1_name = key
		if file_2.endswith(key):
			file_2_name = key

	print(file_1_name)
	print(file_2_name)
	print(corr_value)

	fig , plot = plt.subplots()
	plot.scatter(corr_1, corr_2, s = 5)
	plt.xlabel(file_1_name, fontsize = 20)
	plt.ylabel(file_2_name, fontsize = 20)
	plt.show()

def getCombinations(files):
	combination = []
	for x, y in itertools.combinations(range(len(files)), 2):
		temp = [x,y]
		combination.append(temp)
	return combination


special_case = specialCase(file_name)

samples = getSamples(file_name[0])

for i in range(len(file_name)):
	files[i] = formatDf(file_name[i])

combinations = getCombinations(files)

for pair in combinations:
	if pair[0] == special_case or pair[1] == special_case:
		PCCrow(files[pair[0]], files[pair[1]], file_name[pair[0]],file_name[pair[1]])
	else:
		PCCcol(files[pair[0]], files[pair[1]], samples, file_name[pair[0]], file_name[pair[1]])

