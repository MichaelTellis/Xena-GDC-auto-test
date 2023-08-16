'''
This script checks the accuracy of Xena ETL genomic data from the GDC by independently importing genomic data from the GDC 
and comparing it to the Xena ETL data. This is done by reading the samples of the Xena file and then sending a request to 
the GDC for the file names of the files of each sample. File names which include genomic data are saved into a list. These
file names are then included in another request to the GDC which retrieves each of the file ids. Genomic data of these samples
is imported from the GDC by using these file ids. The data is downloaded as a bundle which are unzipped into various TSV files
each with all genomic data types of a single sample. The raw data in each file is put into a dataframe and formatted based on
Xena requirements. Based upon the data type selected through the arguments, the script compares each cell of the GDC imported 
data and the Xena data to make sure they are identical. If all samples pass 100%, the script will print out success, if not 
the process will stop at where the error is. 



Arguments: 

xena file (arg[1]): 
Use xena gene expression format. Give the asolute path of file. 

available datatypes (arg[2]): 
fpkm_unstranded,
fpkm_uq_unstranded,
tpm_unstranded,
unstranded     (This is star counts)

debug mode (arg[3]): 
When debug mode is True, sample names, file names, ids, will all be printed out. 

'''

import requests
import json
import pandas as pd
import numpy as np
import subprocess
import math
from os import path
import os
import re
import tarfile
import sys

if len(sys.argv[:]) != 4:
    print("Args: this_file   xena_file_path   data_type   debug_mode(True or False)")
    sys.exit(0)
xena_file = sys.argv[1]
if sys.argv[2] == "fpkm":
    data_type = "fpkm_unstranded"
elif sys.argv[2] == "fpkm_uq":
    data_type = "fpkm_uq_unstranded"
elif sys.argv[2] == "tpm":
    data_type = "tpm_unstranded"
elif sys.argv[2] == "star_counts":
    data_type = "unstranded"
else:
    print("Available data types are 'fpkm', 'fpkm_uq', 'tpm', 'star_counts'")
    sys.exit(0)
debug = False
if sys.argv[3] == "True" or sys.argv[3] == "true":
    debug = True


'''
getSamples reads the first line of the xena tsv file, and retrieves all samples into a list (samples), in order.  
'''

def getSamples():
    file = open(xena_file, "r+")
    l = file.readline()
    file.close()
    l = l.split("\t")
    l.pop(0)
    sampleList = []
    for i in l:
        sampleList.append(i.strip())
 
    return(sampleList)
'''
unpeels unneccessary dictionaries and lists returned by GDC API. resJson is the GDC API response in Json format. However,
after being unpeeled it is not in Json format. 
'''

def unpeel(resJson):
    resJson = json.loads(resJson)
    resJson = resJson.get("data")
    resJson = resJson.get("hits")
    return resJson

'''
findFile uses the cases endpoint to retrieve the file names of each sample through a post request. The response is unpeeled. 
Each sample has its own dictionary of file names as well as id. First we search , so for each file name we search for
rna_seq.augmented_star_gene_counts.tsv'. If a file name contains this string, then it will be appended to a new list called 
file_name. 
'''
def findFile(samples): # samples is the list of samples from the xena_file
    cases_endpt = "https://api.gdc.cancer.gov/cases" # we use the cases endpoint for this search and retrieval.
    fields = "files.file_name" # The field we are searching for is the file name
    search = 'rna_seq.augmented_star_gene_counts.tsv' # search is a string that we check in each file name we get, and if it matches, the file name is put in a list called file_list.
    file_list = [] 
    keys = ["files"] # We only want the file names from our response and so we only search inside values with keys named 'files'

    filters = {
        "op": "in",
        "content":{
            "field": "samples.submitter_id",
            "value": samples
            }
    }

    params = {
        "filters": json.dumps(filters),
        "fields": fields,
        "format": "json",
        "size": "20000"
    }
    
    response = requests.post(cases_endpt, headers = {"Content-Type": "application/json"},  json = params)

    #response is the Json response we get from the GDC after sending a POST request.

    # responseJson is our response in a more readable Json format. Once unpeeled however, it is a list containing
    # the id and file names of each sample.
    responseJson = json.dumps(response.json(),indent=2)
    if debug == True:
        print("\n Json formatted response of all file names: ")
        print(responseJson)
    responseJson = unpeel(responseJson)
    
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
    return file_list

'''
getFile is a function that uses the list of file names retrived from the findFile function and sends another request to the GDC
for their file ids. 
'''    

def getFile(file): # file is a list of file_names that we want ids of

    file_endpt = "https://api.gdc.cancer.gov/files" #We use the files endpoint this time

    filters = {
        "op": "in",
        "content":{
            "field": "file_name",
            "value": file
            }
    }

    fields = "file_id"
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
    if debug == True:
        print("\n Json formatted response of all file ids: ")
        print(responseJson)
    responseJson = unpeel(responseJson)

    file_id_list = [] # file_id_list will be a list containing all the file ids.
    for i in responseJson: # i is a dictionary containing a key called 'file_id' and a value that is the file id. 
        file = i["file_id"]
        file_id_list.append(file) # we append the value of i to file_id_list. 

    return file_id_list


# XenaFormat converts the already formatted xena file provided by the user into a xena formatted data frame.

def xenaFormat(xena_file): 
    xena_format = pd.read_csv(xena_file ,sep = '\t', index_col=0)
    xena_format = xena_format.round(10)
    return xena_format


#compare files downloads all of the gene expression files in the selected samples by using the file's id. After extracting them all,
#compare files then goes through each sample's data and compare their values to the provided xena data. If every cell matches, it
#will move on to the next sample until complete. If all values match the program considers it a success, however if there is a difference
# the loop will break.
def compareFiles():
    # Payload is list of file ids in Json format 
    payload = {"ids": file_uuid}

    if debug == True: 
        print("\n file_uuids:")
        print(file_uuid)

    with open("request.txt", "w") as request: # writes payload to a txt file called request.txt
        request.write(str(payload).replace("\'", "\"")) # however, all single quotations must be replaced with double quotations for the download to work. 
        if debug == True:
            print("\n Payload format: ")
            print(request)
    # request.txt includes all of file ids that will be used in a download request of those files. 
    #They are to be downloaded to the file: gdc_download.tar.gz
    #once extracted, this file is a directory with nested directories named after file ids of each file. Inside these nested directories 
    # is one file for each directory. 
    print("\n Importing from GDC: ")
    subprocess.run(["curl", "-o", "gdc_download.tar.gz", "--remote-name", "--remote-header-name", "--request", "POST", "--header", "Content-Type: application/json", "--data", "@request.txt", 'https://api.gdc.cancer.gov/data'])  
    gdc_download = tarfile.open("gdc_download.tar.gz")

    gdc_download.extractall("gdc_download")
    gdc_download.close()
    marker = 0 
    samples_passed = 0 # samples passed increases by one for every sample (file) checked.
    for uuid in file_uuid:

        compare_values = 0  # compare values keeps track of the number of cells compare and if it matches the length of the dataframe
        sample_name = samples[marker]
        # format data downloaded and convert into pd dataframe. First need to find the file using its path.
        file = os.listdir("gdc_download/" + uuid)
        data = pd.read_csv("gdc_download/" + uuid + "/" + file[0], sep ='\t', skiprows = 1, index_col=0)
        data[data_type] = data[data_type] + 1
        gdc_data = pd.DataFrame(np.log2(data[data_type]))
        gdc_data = gdc_data.round(10)
        gdc_data = gdc_data.drop(['N_unmapped', 'N_multimapping', 'N_noFeature', 'N_ambiguous'])
        if debug == True:
            print("\n sample name:")
            print(sample_name)
            print("\n file name: ")
            print(file)
            print("\n file id: ")
            print(uuid)
            print("\n imported data from the GDC dataframe of specified datatype: ")
            print(gdc_data)
        # Compare cells in GDC data with Xena data , if any do not match, loop will break.
        for item in range(len(gdc_data)):
            cell1 = gdc_data.iloc[item][data_type]
            cell2 = xena_data.iloc[item][sample_name]
            if cell1 == cell2:
                compare_values = compare_values + 1
            elif pd.isna(gdc_data.iloc[item][data_type]) and pd.isna(xena_data.iloc[item][sample_name]):
                compare_values = compare_values + 1
            else:
                print("Sample Id: ")
                print(sample_name)
                print("\n GDC download value: ")
                print(gdc_data.iloc[item][data_type])
                print("\n gene id:  ")
                print(gdc_data.loc[item]['gene_id'])
                print("\n Xena format value:" )
                print(xena_data.iloc[item][sample_name])
                print("\n Ensembl ID: ")
                print(xena_data.loc[item]['Ensembl_ID'])
                print("\n Row number: ")
                print(item)
                break
                break
        marker = marker + 1
        if compare_values == len(gdc_data):
            samples_passed = samples_passed + 1
            if debug == True:
                print("\n cells compared: ")
                print(compare_values)
        else:
            print("\nerror or data does not match")
    return samples_passed

samples = getSamples() # samples is a list of all samples in the Xena file
if debug == True:
    print("\n arguments: ")
    print(sys.argv[:])
    print("\n list of all samples: ")
    print(samples)

file_name = findFile(samples) #file_name is a list of all the names of the files we want to download. 
if debug == True:
    print("\n list of all file names (they are not in order): ")
    print(file_name)
file_uuid = getFile(file_name) # file_uuid is the list of all of the ids of the files we want to download. 
if debug == True:
    print("\n list of all file ids: ")
    print(file_uuid)
xena_data = xenaFormat(xena_file) # xena_data is a data frame of the provided data.
if debug == True:
    print("\n \n Xena format DataFrame: \n")
    print(xena_data)
samples_pass = compareFiles() # samples_pass is the number of samples that were 100% accurate.

if samples_pass == len(samples): # checks if the number of smaples passed is equal to the number of samples. 
    print("\n samples compared: ")
    print(samples_pass)
    print("\nSuccess")