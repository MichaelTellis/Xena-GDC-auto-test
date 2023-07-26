import requests
import json
import pandas as pd
import numpy as np
import subprocess
import math
from os import path
templateFile =  "/Users/michaeltellis/XenaGDC/CTSP-DLBCL1.star_fpkm-uq.tsv"

def getFields():
    file = open(templateFile, "r+")
    l = file.readline()
    file.close()
    l = l.split("\t")
    l.pop(0)
    sampleList = []
    for i in l:
        sampleList.append(i.strip())
 
    return(sampleList)

def unpeel(resJson):
    resJson = json.loads(resJson) # transforms Json into dictionary
    # unpeels nested dictionaries and lists
    resJson = resJson.get("data") 
    resJson = resJson.get("hits")
  #  resJson = resJson[0]
    return resJson


def findFile(samples):
    cases_endpt = "https://api.gdc.cancer.gov/cases"

    filters = {
            "op": "in",
            "content":{
                "field": "samples.submitter_id",
                "value": samples
                }
    }

    # With a GET request, the filters parameter needs to be converted
    # from a dictionary to JSON-formatted string
    #fields = getFields()
    fields = "files.file_name"
    params = {
        "filters": json.dumps(filters),
        "fields": fields,
        "format": "json",
        "size": "100"
        }
    response = requests.post(cases_endpt, headers = {"Content-Type": "application/json"},  json = params)

    #print(json.dumps(response.json(),indent=2))
    #responseJson has the list of files of the given sample
    responseJson = json.dumps(response.json(),indent=2)


    responseJson = unpeel(responseJson)
    search = 'rna_seq.augmented_star_gene_counts.tsv'
    file_list = []
    keys = ["files"]
    for i in responseJson:
      
        item = []
        for key, value in i.items():
            if key in keys:
                item = value
                for i in item:
                    file = i["file_name"]
                    if file.find(search) != -1:
                        file_list.append(file)
                        break
                break

        
        
        
    return file_list
samples = getFields()
file_name = findFile(samples)
print(file_name)



def getFile(file):

    cases_endpt = "https://api.gdc.cancer.gov/files"

    filters = {
        "op": "in",
        "content":{
            "field": "file_name",
            "value": file
            }
    }


    #fields = getFields()
    fields = "file_id"
    params = {
        "filters": json.dumps(filters),
        "fields": fields,
        "format": "json",
        "size": "100"
        }

    response = requests.post(cases_endpt, headers = {"Content-Type": "application/json"},  json = params)

    #print(json.dumps(response.json(),indent=2))
    #responseJson has the list of files of the given sample
    responseJson = json.dumps(response.json(),indent=2)
    responseJson = unpeel(responseJson)
    search = 'file_id'
    file_id = ""
    for key, item in responseJson.items():
        if search == key:
            file_id = item
            break 
    return file_id
#Download data using GDC api, only downloads once per file.
file_uuid = getFile(file_name)
print("\n file_uuid:")
print(file_uuid)
file_url = "https://api.gdc.cancer.gov/data/" + file_uuid
print("\n Importing from GDC: ")
subprocess.run(["curl", "--remote-name", "--remote-header-name", file_url])

# format data downloaded and convert into pd dataframe
data = pd.read_csv(file_name, sep ='\t', skiprows = 1, index_col=0)

data['fpkm_uq_unstranded'] = data['fpkm_uq_unstranded'] + 1
gdc_star_fpkm_UQ_data = pd.DataFrame(np.log2(data['fpkm_uq_unstranded']))
print("\n \n Imported from GDC DataFrame: ")
print(gdc_star_fpkm_UQ_data)

# convert Xena formated data into pd dataframe
star_fpkm_UQ_data = pd.read_csv(templateFile ,sep = '\t', index_col=0)
star_fpkm_UQ_data = star_fpkm_UQ_data.round(10)
gdc_star_fpkm_UQ_data = gdc_star_fpkm_UQ_data.round(10)
print("\n \n Xena format DataFrame: \n")
print(star_fpkm_UQ_data)

# Compare cells in Downloaded and formated GDC data with Xena data retrieved, if any do not match, loop will break.
compare_values = 0

for item in range(len(gdc_star_fpkm_UQ_data)):
    cell1 = gdc_star_fpkm_UQ_data.iloc[item]['fpkm_uq_unstranded']
    cell2 = star_fpkm_UQ_data.iloc[item]['DLBCL11296-sample']
    if cell1 == cell2:
        compare_values = compare_values + 1
    elif pd.isna(gdc_star_fpkm_UQ_data.iloc[item]['fpkm_uq_unstranded']) and pd.isna(star_fpkm_UQ_data.iloc[item]['DLBCL11296-sample']):
        compare_values = compare_values + 1
    else:
        print("star fpkm-UQ Unstranded GDC format value: ")
        print(gdc_star_fpkm_UQ_data.iloc[item]['fpkm_uq_unstranded'])
        print("\n gene id:  ")
        print(gdc_star_fpkm_UQ_data.iloc[item]['gene_id'])
        print("\n star fpkm-UQ Unstranded Xena format value:" )
        print(star_fpkm_UQ_data.iloc[item]['DLBCL11296-sample'])
        print("\n Ensembl ID: ")
        print(star_fpkm_UQ_data.iloc[item]['Ensembl_ID'])
        print("\n Row number: ")
        print(item)
        break
if compare_values == len(gdc_star_fpkm_UQ_data):
    print("\nSuccess")
else:
    print("\nerror or data does not match")
#compare_values = gdc_star_fpkm_UQ_data.iloc[]['fpkm_uq_unstranded'].isin(star_fpkm_UQ_data['DLBCL11296-sample']).value_counts()
#print(compare_values)
#num_Nan = gdc_star_fpkm_UQ_data['fpkm_uq_unstranded'].value_counts()['']
#print(num_Nan)