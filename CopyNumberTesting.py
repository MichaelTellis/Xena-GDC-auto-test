import pandas as pd
import numpy as np
import json 
import requests
import subprocess
import math
from os import path
import os
import re
import tarfile
import sys


'''
Usage Instructions: 
Command Line:
[This file name] [Xena-GDC-ETL imported .tsv file name]

example: python3 XenaGDC/CopyNumberTesting.py /Users/Downloads/CGCI-HTMCP-LC.cnv_ascat-ngs.tsv
'''



#### CONSTANTS ########################################
# example file:
#xena_file = "/Users/michaeltellis/Downloads/CGCI-HTMCP-LC.cnv_ascat-ngs.tsv"
xena_file = sys.argv[1]
keyword = ["files", "samples", "Chrom", "Start", "End", "value"]
cases_endpt = "https://api.gdc.cancer.gov/cases"
file_endpt = "https://api.gdc.cancer.gov/files"
fileRequestFields = "files.file_name"# The field we are searching for is the file name
fileRequestSearch = 'copy_number_variation.seg.txt' # search is a string that we check in each file name we get, and if it matches, the file name is put in a list called file_list.
fileRequestKeys = ["files"]  # We only want the file names from our response and so we only search inside values with keys named 'files'
downloadFilesFields = "file_id"
orderedSamplesFields = ["samples.submitter_id", "files.file_id"] # The field we are searching for is the file name

### Globals #######################################################################
xena_df = pd.read_csv(xena_file, sep="\t")
sample = []
file_names = []

print(xena_df)
def sample(df):
	samples = df["sample"].tolist()
	sample_list = []
	[sample_list.append(x) for x in samples if x not in sample_list]
	print(sample_list)
	return sample_list


def unpeel(resJson):
    resJson = json.loads(resJson)
    resJson = resJson.get("data")
    resJson = resJson.get("hits")
    return resJson


def fileRequest(fields, search, keys, sample_list, cases_endpt):
	file_list = [] 
	filters = {
		"op": "in",
		"content":{
			"field": "samples.submitter_id",
			"value": sample_list
		}
	}

	params = {
		"filters": json.dumps(filters),
		"fields": fields,
		"format": "json",
		"size": "20000"
	}
	    
	response = requests.post(cases_endpt, headers = {"Content-Type": "application/json"},  json = params)
	responseJson = json.dumps(response.json(),indent=2)
	responseJson = unpeel(responseJson)
	print(responseJson)

	for i in responseJson: # i is a dictionary in the responseJson list which correspondes to a single sample.  
		item = [] 
		for key, value in i.items(): 
			if key in keys: # if this specific key in i is equal to 'files', then its value will be the list, item. 
				item = value # item now is a list of the file names of all files in a specific sample (as well as the sample ids.)
				for i in item:
					file = i["file_name"] # each element in the item list is a dictionary and if that dictionary's key is 'file_name' 
					if file.find(search) != -1: # then we check if the value( file name) matches our search string.
						file_list.append(file) # if it does, it is appended to a new list of the files we want to download. 
						break
				break
	print(file_list)
	return file_list

def downloadFiles(fields, file_endpt, file_list):
	file_id_list = [] # file_id_list will be a list containing all the file ids.
	filters = {
		"op": "in",
		"content":{
			"field": "file_name",
			"value": file_list
			}
	}

	params = {
		"filters": json.dumps(filters),
		"fields": fields,
		"format": "json",
		"size": "100"
		}

	response = requests.post(file_endpt, headers = {"Content-Type": "application/json"},  json = params)
	#response is the Json response we get from the GDC after sending a POST request.

	#response Json is response in a more readable format. 
	responseJson = json.dumps(response.json(),indent=2)
	responseJson = unpeel(responseJson)
	for i in responseJson: # i is a dictionary containing a key called 'file_id' and a value that is the file id. 
		file = i["file_id"]
		file_id_list.append(file) # we append the value of i to file_id_list. 
	print(file_id_list)

	payload = {"ids": file_id_list}

	with open("request.txt", "w") as request: # writes payload to a txt file called request.txt
		request.write(str(payload).replace("\'", "\""))
	subprocess.run(["curl", "-o", "gdc_download.tar.gz", "--remote-name", "--remote-header-name", "--request", "POST", "--header", "Content-Type: application/json", "--data", "@request.txt", 'https://api.gdc.cancer.gov/data'])  
	gdc_download = tarfile.open("gdc_download.tar.gz")

	gdc_download.extractall("gdc_download")
	gdc_download.close()

	return file_id_list



def orderedSamples(cases_endpt, fields,file_id_list):
	fields = ",".join(fields)
	file_list = [] 
	filters = {
		"op": "in",
		"content":{
			"field": "files.file_id",
			"value": file_id_list
		}
	}
		
	params = {
		"filters": json.dumps(filters),
		"fields": fields,
		"format": "json",
		"size": "20000"
	}
		    
	response = requests.post(cases_endpt, headers = {"Content-Type": "application/json"},  json = params)
	responseJson = json.dumps(response.json(),indent=2)
	responseJson = unpeel(responseJson)
	return responseJson

def compareSamples(responseJson, file_id_list, xena_df, sample_list):
	samples_passed = 0
	for i in responseJson:
		marker = None
		idIndex = None
		sample = None
		chromCompare = 0
		startCompare = 0
		endCompare = 0
		valueCompare = 0
		for key, value in i.items():
			if key == keyword[0]:
				for j in value:
					#print(j)
					for name, item in j.items():
						#print(item)
						if item in file_id_list:
							idIndex = file_id_list.index(item)
							break
							break

		for key, value in i.items():
			if key == keyword[1]:
				for j in value:
					for name, item in j.items():
						if item in sample_list:
							sample = item
							break 
							break
		uuid = file_id_list[idIndex]
		file = os.listdir("gdc_download/" + uuid)
		data = pd.read_csv("gdc_download/" + uuid + "/" + file[0], sep ='\t', skiprows = 1, index_col=0, header=None)
		data.columns = ["Chrom", "Start", "End", "value", "", ""]
		pos = []
		result = xena_df.isin([sample])
		seriesObj = result.any()
		columnNames = list(seriesObj[seriesObj == True].index)

		for col in columnNames:
			rows = list(result[col][result[col] == True].index)
			for row in rows:
				pos = [row, col]
				break
			break
		#print(pos)

		
		zeroCell = xena_df.loc[pos[0]][pos[1]]
		print(zeroCell)
		for item in range(len(data)):
			cell2Pos = pos[0] + item
			cell1 = data.iloc[item][keyword[2]]
			cell2 = xena_df.loc[cell2Pos][keyword[2]]
			if cell1 == cell2:
				chromCompare += 1
			else:
				print("fail")
				print(cell1)	
				print(cell2)
				print(keyword[2])
				break
		for item in range(len(data)):
			cell2Pos = pos[0] + item
			cell1 = data.iloc[item][keyword[3]]
			cell2 = xena_df.loc[cell2Pos][keyword[3]]
			if cell1 == cell2:
				startCompare += 1
			else:
				print("fail")
				print(cell1)	
				print(cell2)
				print(keyword[3])
				break
		for item in range(len(data)):
			cell2Pos = pos[0] + item
			cell1 = data.iloc[item][keyword[4]]
			cell2 = xena_df.loc[cell2Pos][keyword[4]]
			if cell1 == cell2:
				endCompare += 1
			else:
				print("fail")
				print(cell1)	
				print(cell2)
				print(keyword[4])
				break
		for item in range(len(data)):
			cell2Pos = pos[0] + item
			cell1 = data.iloc[item][keyword[5]]
			cell2 = xena_df.loc[cell2Pos][keyword[5]]
			if cell1 == cell2:
				valueCompare += 1
			else:
				print("fail")
				print(cell1)	
				print(cell2)
				print(keyword[5])
				break
		if chromCompare == len(data) and startCompare == len(data) and endCompare == len(data) and valueCompare == len(data):
			samples_passed += 1
			print("pass")
		else: 
			print("fail")

		print(chromCompare)
		print(startCompare)
		print(endCompare)
		print(valueCompare)
	if samples_passed == len(sample_list):
		print("success")
	else: 
		print("fail")



sample = sample(xena_df)

file_names = fileRequest(fileRequestFields, fileRequestSearch, fileRequestKeys, sample, cases_endpt)

file_uuids = downloadFiles(downloadFilesFields, file_endpt, file_names)

orderedSamples = orderedSamples(cases_endpt, orderedSamplesFields, file_uuids)

compareSamples(orderedSamples, file_uuids, xena_df, sample)

