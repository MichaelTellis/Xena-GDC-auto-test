import requests
import json
import sys
import pandas as pd
import numpy as np
import sys
'''
Xena Survival Analysis Endpt Validation works very similarily to Xena Survival Matrix Validation. The major difference 
between the two is that this script uses the /analysis/survival endpoint instead of the /cases endpoint to retrieve
survival data. First, a post request is sent to the endpoint for all survival data of a specific project. Interestingly, 
the fields entered in the post request does not affect the response, which is why 'fields' is empty. Then, we loop 
through the Json response and add any OS_time data to 'time' list and submitter ids to 'submitter_id'. The json response
also includes OS data. However, it was decided to import the OS data independently using the cases endpoint. However, 
the json response recieved is not in the same order as the list of submitter ids sent, so all OS data has to be matched 
with other survival data. Once this is done, all survival data is formatted into a dataframe and the same is done
for the xena TSV file. The two dataframes are then compared for any missing data or incorrect data.

To run this script, 3 arguments are required: this file's name and path, the xena file's name and
path, and the project id. 
'''
############################ Globals #####################################

time = [] # Where all values for the of events will be tracked. (OS_time)

status = []  # Where all values for the vital_status of cases will be tracked (OS)

submitter_id = [] # The submitter ids of all cases that will be compared

gdc_df = pd.DataFrame() # data frame of all survival data directly imported from the GDC.

xena_df = pd.DataFrame() # data frame of all survival data in the Xena file given.


########################## Constants #####################################


if len(sys.argv) != 3:
	print("incorrect arguments")
	sys.exit(0)

xena_file = sys.argv[1] # Xena file path and name, second argument

project_id = sys.argv[2] # project id, third argument. 

#list of keywords that the getData() function will search for in order to correctly assign values to the time, status, or submitter_id lists. 
survival_keys = ["time", "submitter_id"]

#/analysis/survival endpoint link
survival_endpt = "https://api.gdc.cancer.gov/analysis/survival"

#/cases endpoint link
cases_endpt = "https://api.gdc.cancer.gov/cases"

# fields that are requested by the getStatus() function
status_fields = ["demographic.vital_status", "submitter_id"]

# filter used by the getStatus() function
status_filter = "submitter_id"
# list of keywords that the getStatus() function will search in the json response in order to find all submitter ids.
status_keys = ["submitter_id","demographic", "vital_status"]

# The column names of the 3 columns in 'gdc_df'.
column_names = ["OS.time", "OS", "_PATIENT"]


'''
unpeel is a function that will unpeel the nested dictionaries that are in the GDC json response. 
'''
def unpeel(resJson):
	resJson = json.loads(resJson) # resjson is a json formatted dictionary
	resJson = resJson.get("results")
	resJson = resJson[0]
	resJson = resJson.get("donors")
	return resJson # however, it is returned as a list

'''
getData sends a post request to the GDC /analysis/survival endpoint for all survival data of a specific project. 
In this case fields are not needed as the response alreay includes all required data. After the json is enpeeled,
a for loop iterates through each element (case) in the response json. If the key matches the 'time' keyword, the 
value is appended to the 'time' list. If it matches the 'submitter_id' keyword it is appended to the 'submitter_id'
list.
'''
def getData(endpt, keys, project): # /analysis/survival endpoint is used. Project is project id. keys are the keywords 
	# searched in responseJson
	fields = []
	fields = ",".join(fields)
	filters = {
        "op": "in",
        "content":{
            "field": "cases.project.project_id",
            "value": project
			}
        
	}

	# A POST is used, so the filter parameters can be passed directly as a Dict object.
	params = {
    	"filters": filters,
    	"fields": fields,
    	"format": "json",
    	"size": "2000"
    	}

	# The parameters are passed to 'json' rather than 'params' in this case
	response = requests.post(endpt, headers = {"Content-Type": "application/json"}, json = params)

	responseJson = json.dumps(response.json(),indent=2) #Json formatted
	responseJson = unpeel(responseJson) #unpeeled into a list
	for i in responseJson: # iterate through element inside list (each is a case)
		for key, value in i.items(): # each element is a dictionary that we iterate through as well
			if keys[0] in key: # if the key matches the keyword: 'time', the value is appended to time
				time.append(value)
			elif keys[1] in key:  # if the key matches the keyword: 'submitter_id', the value is appended to submitter_id
				submitter_id.append(value)

	print(time)
	print(submitter_id)

'''
getStatus uses the cases endpoint in order to independently retrieve OS data from all cases of a project. After getting the json 
response and turning it into a list, we iterate through each element in the list. First we find the submitter_ids that we will use 
for mapping the data to match data from getData(). Then we iterate again in order to get the OS data, which is put in the vital_status
list in the correct order.
'''
def getStatus(endpt, fields, submitter_id, filter_field, keys): #cases endpt, fields: ["demographic.vital_status", "submitter_id"]
	#'submitter_id' is list of submitter ids from getData, filter_field: 'submitter_id' is the field we're filtering by, keys:
	#key words searched in the response json: ["submitter_id","demographic", "vital_status"]

	vital_status = [None] * len(submitter_id) #vital_status becomes a list the length of submitter_id where each element is None.
	submitter_id_dict = dict.fromkeys(submitter_id) # turns the submitter_id list into dict format. 
	index = 0 # will be used to put an index as the value
	for key, value in submitter_id_dict.items():
		submitter_id_dict[key] = index # each value is given an index which increases by one as we go through the loop.
		index +=1
	fields = ",".join(fields) #fields need to be a comma seperated string
	filters = {
        "op": "in",
        "content":{
            "field": filter_field, 
            "value": submitter_id
			}
        
	}

	# A POST is used, so the filter parameters can be passed directly as a Dict object.
	params = {
    	"filters": filters,
    	"fields": fields,
    	"format": "json",
    	"size": "2000"
    	}

	# The parameters are passed to 'json' rather than 'params' in this case
	response = requests.post(endpt, headers = {"Content-Type": "application/json"}, json = params)

	responseJson = json.dumps(response.json(),indent=2)
	#responseJson is a json format dict. 
	responseJson = json.loads(responseJson)
	responseJson = responseJson.get("data")
	responseJson = responseJson.get("hits")
	#responseJson is now a list
	for i in responseJson: # iterate through each element (case) in list. Each elements is a dictionary.
		mapping_index = None # create mapping_index variable
		for key, value in i.items(): #iterate through said dictionary
			if keys[0] in key: # checks if keyword matches dictionary key, in this case 'submitter_id'. 
				if value in submitter_id_dict.keys(): #Makes sure that the submitter id matches a submitter_id in 'submitter_id_dict'.
					mapping_index = submitter_id_dict[value] #If so, find the value at said submitter id. This value is an index that tells us which position the OS should be.
					#mapping_index will contain the correct position this OS should be. 
					break
		for key, value in i.items(): #iterate through the same dictionary a second time
			if keys[1] in key: # checks if keyword matches dict key, in this case 'demgraphic'				
				for name, item in value.items(): # unlike 'submitter_id', 'demographic' has a dict within it. which is now iterated through.
					if keys[2] in name and item == "Dead":  #If dict key matches 'vital_status' and if the value is 'Dead'. 
						vital_status[mapping_index] = 1		# if so, then we know the OS value should be one, and we use mapping_index to get correct position
					elif keys[2] in name and item == "Alive": #If dict key matches 'vital_status' and if the value is 'Alive'.  OS must be 0.
						vital_status[mapping_index] = 0
					else: 
						print("error in vital status") 
						break
						break
						break
	print(vital_status)
	return vital_status

'''
Now all required data is collected, it is formatted into a dataframe, where the first column is OS_time, second is OS, third is submitter_id
The dataframe is soted by OS_time.
'''
def formatData(col):
	gdc_data = {col[0]:time, col[1]: status, col[2]: submitter_id} # col is a list of column names, found in constants as 'column_names'
	df = pd.DataFrame(data=gdc_data)
	df = df.sort_values(by = col[0], ascending = True)
	return df #returns to gdc_df

'''
The same is done with the xena TSV file 
'''
def xenaFormat(file_name):
	xena_format = pd.read_csv(file_name, sep= '\t')
	return xena_format # xena_format is dataframe
'''
comparison compares both dataframes. The GDC imported dataframe is named 'gdc', the xena dataframe is named 'xena'. However, the two dataframes are in different orders.
So the first thing that is done is to match each row (case) in the GDC data to a row in the Xena data. This is done by using the submitter_ids to map one case onto another. Afterwards,
the two dataframes are compared on a cell to cell basis. If a cell in the 'OS_time' column matches, a marker for OS_time increases by one. The same is true for the 'OS' column. Once all 
cells have been compared, the 'OS_time_compare' is compared with 'OS_compare' and the length of 'gdc'. If all three are equal, the program is successful, otherwise, the program fails. 
A list of missing cases, their submitter_id, OS_time, and OS is also printed. 
'''
def comparison(gdc, xena, column):
	os_time_compare = 0 #used to track successful OS_time comparisons
	os_compare = 0 #used to track successful OS comparisons.
	missing_list = [] # cases that are missing from the xena file are appended to this list.
	for cell in range(len(gdc.index)): # cell is a marker that goes through each row of 'gdc'. 
		search = gdc.loc[cell, column[2]] # search is a string that corresponds to specified row's submitter_id
		if (xena[column[2]].eq(search)).any() == True: # check if that submitter_id exists in the 'xena' dataframe. If so, compare the other two columns. 
			row, col = np.where(xena == search) #finds row and column where submitter_id is located in 'xena'
			xena_time = xena.loc[row[0], column[0]] #xena_time corresponds to the OS_time at specified row and column named 'OS_time'
			gdc_time = gdc.loc[cell, column[0]] # gdc_time corresponds to the OS_time value at row 'cell'.
			xena_status = xena.loc[row[0], column[1]] # xena_status is the OS value at the specified row (based on submitter_id)
			gdc_status = gdc.loc[cell,column[1]]# gdc_status corresponds to the OS value at row 'cell'.
			if xena_time == gdc_time: #check if the two OS_times are equivalant. 
				os_time_compare += 1
			if xena_status == gdc_status: # check if the two OSs are equivalent.
				os_compare += 1
		else:   # if 'search' (submitter_id) does not exist, then the case is missing from the xena file
			print("\nCase is missing from survival data matrix (submitter_id):")
			print(search)
			print("\n" + column[0] + " value:")
			print(gdc.loc[cell, column[0]])
			print("\n" + column[1] + " value:")
			print(gdc.loc[cell, column[1]])
			case_info = [search, gdc.loc[cell, column[0]], gdc.loc[cell, column[1]]] # the survival data on the case is then saved to 'missing_list'.
			missing_list.append(case_info)
	if os_compare == len(gdc.index) and os_time_compare: # Checks if OS_time_compare and OS_compare are equivalent to the length of the 'gdc' dataframe.
		print("\nsuccess!")
		print("number of cases successfully compared: " + str(os_compare))
	else:     # if they are not equivalent the program fails. 
		print("\nFailed")
		print("\nCases successfully compared:" + str(os_compare) + "/" + str(len(gdc.index)))
		print("\nList of missing cases [Submitter_id, OS_time, OS]:")
		print(missing_list)
		print("\nNumber of cases where vital status was identitcal: " + str(os_compare))
		print("\n number of cases where time was identitcal: " + str(os_time_compare))

#gets the OS_time data and submitter_ids of a project using a post request.
getData(survival_endpt, survival_keys, project_id)

# gets the OS data and matches it in the same order as getData()
#stats is the ordered list of OS data.
status = getStatus(cases_endpt, status_fields, submitter_id, status_filter, status_keys)

# all of this data is then formatted into a Dataframe called 'gdc_df'
gdc_df = formatData(column_names)
# the xena file TSV is also formatted into a dataframe called 'xena_df'
xena_df = xenaFormat(xena_file)
print(xena_df)
print(gdc_df)

# The two dataframes are compared and any missing cases or incorrect data is printed out at the end. 
comparison(gdc_df, xena_df, column_names)

