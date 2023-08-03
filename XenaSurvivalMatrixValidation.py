import requests
import json
import sys
import pandas as pd
import numpy as np
'''
Xena survival matrix validation is a script that checks the accuracy of data imported from the GDC 
(Genomic Data Commons) by the Xena ETL code. It uses an independent method to retrieve survival data
from a specified project in the GDC. While the Xena ETL code uses the /analysis/survival endpoint
to retrieve surival data, this program uses the /cases endpoint. The /cases endpoint has multiple 
fields relevant to survival data (found in the list survival_fields under the constants section).
After sending a post request to the GDC API that asks for the submitter ids of the selecred project,
the script recieves a json response. We then loop through the json to get the list of all submitter ids. 
A second post request is created with the previously retrieved submitter ids as filters and the
survival_fields list as fields. We then loop through the recieved json and assign submitter ids, 
event times, and vital status to different lists. These lists will become 3 columns in a data 
frame named 'survival_df'. However, if the survival time of a case is less than or equal to 0,
the case is 'thrown out' as it will not be useful when analyzing such data. A second data frame
named 'xena_df' is created by reading the provided xena formatted tsv file to the data frame. 
Finally, by using the submitter ids as indexes, the two data frames are compared for any missing 
cases in the xena data frame. If there is, the submitter id, OS_time, and OS will be added to a
list named 'missing_list'. Otherwise, the program will return a success, meaning all data matches. 

To run this script, 3 arguments are required: this file's name and path, the xena file's name and
path, and the project id. 
'''

#################Global: ######################################################################

all_time = [] # Where all values for the of events will be tracked.

all_status = [] # Where all values for the vital_status of cases will be tracked

all_submitter_id = [] # The submitter ids of all cases that will be compared

submitter_id_list = [] # A secondary list of submitter ids that is only used for the getSubmitterId() function

survival_df = pd.DataFrame() # data frame of all survival data directly imported from the GDC.

xena_df = pd.DataFrame() # data frame of all survival data in the Xena file given.

################## Constants: ##############################################################

if len(sys.argv) != 3:
    print("incorrect arguments")
    sys.exit(0)

xena_file = sys.argv[1] # Xena file path and name, second argument

project = sys.argv[2] # project id, third argument. 
# fields that are requested by the getData() function
survival_fields = ["demographic.days_to_death", "samples.submitter_id", 
"submitter_id", "diagnoses.days_to_best_overall_response",
"diagnoses.days_to_last_follow_up", 
"diagnoses.days_to_last_known_disease_status", 
"diagnoses.days_to_recurrence",
"follow_ups.days_to_adverse_event", 
"follow_ups.days_to_comorbidity",
"follow_ups.days_to_follow_up",
"follow_ups.days_to_progression",
"follow_ups.days_to_progression_free",
"follow_ups.days_to_recurrence",
"demographic.vital_status"]

#list of keywords that the getData() function will search for in order to correctly assign values to the all_time, all_status, or all_submitter_id lists. 
survival_keys = ["submitter_id", "samples", "demographic", "diagnoses", "follow_ups", "Alive", "Dead"]

# filter used by the getData() function. 
survival_filter = "submitter_id"

# fields that are requested by the getSubmitterId() function
submitter_id_fields = ["submitter_id", "case_id"]

# list of keywords that the getSubmitterId() function will search in the jsonhttps://www.youtube.com/ response in order to find all submitter ids.
submitter_id_keys = ["submitter_id"]

# filter used by the getSubmitterId() function
submitter_id_filters = "project.project_id"

#cases endpoint
cases_endpt = "https://api.gdc.cancer.gov/cases"

# The column names of the 3 columns in 'survival_df'.
column_names = ["OS.time", "OS", "_PATIENT"]



'''
unpeel is a function that will unpeel the nested dictionaries that are in the GDC json response. 
'''
def unpeel(resJson):
    resJson = json.loads(resJson) # resjson is a json formatted dictionary
    resJson = resJson.get("data")
    resJson = resJson.get("hits")
    return resJson # however, it is returned as a list

'''
The first step of the program is to get all of the submitter ids from the project. These are very important as each of
them match a single case in the specified project. Project is the project id, keys is a list of keywords that the for loop 
searches for in the response json. submitter_id_fields currently contains the field 'submitter_id' as we only need the submitter id.
filter_field is just the field we will use to filter. It is currently 'projects.project_id' as we are using the project_id as our filter. 
The first part of the function sets up a post request in json format. The second part of the function are several for loops that go 
through each element of the response list to search for the submitter id and append it to the list 'submitter_id' whic his returned at the 
end of the funciton.
'''

def getSubmitterId(project, keys, submitter_id_fields, filter_field):
    # all submitter ids recieved from the request will be appended to submitter_id 
    submitter_id = []
    # the fields must be a string joined together by commas
    fields = ",".join(submitter_id_fields) 

    filters = {
        "op": "in",
        "content":{
            "field": filter_field,
            "value": project
            }
    }

    # A POST is used, so the filter parameters can be passed directly as a Dict object.
    params = {
        "filters": filters,
        "fields": fields,
        "format": "json", # json format can also be TSV format
        "size": "2000"
        }

    # The parameters are passed to 'json' rather than 'params' in this case
    response = requests.post(cases_endpt, headers = {"Content-Type": "application/json"}, json = params)
    # response json is a readable dictionary in json format
    responseJson = json.dumps(response.json(),indent=2)
    responseJson = unpeel(responseJson)
    #following the unpeel() it is a list where each element is a single case that carries its uuid and submitter id. 

    for i in responseJson:  # response is a list so i is a dictionary that is one element in that list. This dictionary correspondes to one case.
        for key, value in i.items(): # we go through each key-value pair in the dictionary.
            if key in keys: # and search if the key matches 'keys' which is currently "submitter_id"
                submitter_id.append(value) # if so, then the value is the submitter_id and is appended to the submitter_id list. 
                break
    print(submitter_id)
    return(submitter_id)


'''
Once there is a list of submitter_ids. We can use that to find anything we want about each case. In this case, we want to 
find survival data. We pass in the submitter_ids in the list 'submitter_id', keywords to filter through the response json 
are in the list 'keys', fields we want data for are in the list 'survival_fields', and 'filter_field' is a field we will filter by,
(in this case 'submtter_id'). This is then packaged into a post request which returns with a json. Json is then passed into the nested for loops 
in order to split data into three lists, 'all_submitter_id' (submitter_ids), 'all_time' (time events), 'all_status' (vital_status where 1 is dead and 0 is alive).
'''
def getData(submitter_id, keys, survival_fields, filter_field):
    # the fields must be a string joined together by commas
    fields = ",".join(survival_fields)
    
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
        "format": "json", # json format can also be tsv format
        "size": "2000"
        }

    # The parameters are passed to 'json' rather than 'params' in this case
    response = requests.post(cases_endpt, headers = {"Content-Type": "application/json"}, json = params)
    #response json is a readable dictionary in json format
    responseJson = json.dumps(response.json(),indent=2)
    responseJson = unpeel(responseJson)# however once passed through unpeel(), it is a list
    # i is an element in responsejson and is dictionary type. It corresponds to a single case
    for i in responseJson:
        # time, status, and submitter id will include the times, vital_status and submitter_id found in i.
        time = [] 
        status = []
        submitter_id = []
        # we go through each key values pair in the dictionary. 
        for key, value in i.items(): 
            if keys[0] in key: # keys[0] is a keyword, in this case it is 'submitter_id'
                submitter_id.append(value) # if the key is 'submitter_id' then the value is the submitter_id of the case and is appended to the correct list.
            elif key in keys[2]:  # keys[2] is a keyword, in this case it is 'demographic'.
                for name, item in value.items(): # however, demographic has several fields under it which are within another nested dict. 
                    if item != None: # if there is no data on demographics, the value in the dict will be None, if its not, we check it. 
                        if item in keys[5:7]: # if the value is 'Alive' or 'Dead' then we know it is vital status and the value is appended to the correct list. 
                            status.append(item)
                        else: #Otherwise the value is appended to 'time' since all other values are related to OS_time. OS_time field: ('demographic.days_to_death')
                            time.append(item)
                            if item == 0:
                                print(name)
                                print(item)
                                print(submitter_id)
                            
            elif key in keys[3]: # keys[3] is a keyword, in this case it is 'diagnoses'.
                for name, item in value[0].items(): # however, all the data that is relavant is within the first element of diagnoses, which is also a nested dict.
                    if item != None:# if there is no data on diagnoses the value in the dict will be None, if its not, we check it. 
                        if item == keys[5:7]:  # if the value is 'Alive' or 'Dead' then we know it is vital status and the value is appended to the correct list. 
                            status.append(item)
                        else:  #Otherwise the value is appended to 'time' since all other values are related to OS_time. 
                            # OS_time fields: (diagnoses.days_to_best_overall_response", "diagnoses.days_to_last_follow_up", "diagnoses.days_to_last_known_disease_status", 
                            #"diagnoses.days_to_recurrence",)
                            time.append(item)
                            if item == 0:
                                print(name)
                                print(item)
                                print(submitter_id)
                            
            elif key in keys[4]: # keys[4] is a keyword, in this case it is 'follow_ups'
                for j in value: # the value for 'follow_ups' is a list that we must iterate through
                    for name, item in j.items(): # we check the dictionary inside the list
                        if item!= None: # if there is no data on follow_ups the value in the dict will be None, if its not, we check it. 
                            if item == keys[5:7]: # if the value is 'Alive' or 'Dead' then we know it is vital status and the value is appended to the correct list. 
                                status.append(item)
                            else:   #Otherwise the value is appended to 'time' since all other values are related to OS_time. 
                                #OS_time fields: "follow_ups.days_to_adverse_event", "follow_ups.days_to_comorbidity","follow_ups.days_to_follow_up",
                                #"follow_ups.days_to_progression", "follow_ups.days_to_progression_free", "follow_ups.days_to_recurrence"
                                time.append(item) 
                                if item == 0:
                                
                                    print(name)
                                    print(item)
                                    print(submitter_id)
                                
        if len(time) > 1: # sometimes, there are multiple values appended to the 'time' list, so only the Max is taken.
            time = [max(time)]
        if len(time) == 1: # Other times, there is no time (len of 0), which also means no survival data so only those with time data are included.

            if status[0] == "Alive": # changing the OS from a string to a 1 or 0. 
                status[0] = 0
            elif status[0] == "Dead":
                status[0] = 1   
            if time[0] > 0: # times of 0 or less are also not relevant for survival data and are not included. 
                all_time.append(time[0]) 
                all_status.append(status[0])
                all_submitter_id.append(submitter_id[0])

'''
Once all of the data has been corectly filtered, it is formatted into a data frame. The first column has all OS_time values
the second column has all OS values, the third column has all submitter_ids. Each of them also have identical names to the 
Xena column names. The values are then sorted in ascending order based on OS_time. 
'''

def formatData(col):
    
    gdc_data = {col[0]:all_time, col[1]: all_status, col[2]: all_submitter_id}
    gdc_df = pd.DataFrame(data=gdc_data)
    gdc_df = gdc_df.sort_values(by = col[0], ascending = True)
    return gdc_df # Dataframe returned
'''
The Xena file also needs to be formatted to a data frame. Since it is a TSV file, that we use the read_csv function.
'''
def xenaFormat(file_name):
    xena_format = pd.read_csv(file_name, sep= '\t')
    return xena_format # Dataframe returned

'''
comparison compares both dataframes. The GDC imported dataframe is named 'survival', the xena dataframe is named 'xena'. However, the two dataframes are in different orders.
So the first thing that is done is to match each row (case) in the GDC data to a row in the Xena data. This is done by using the submitter_ids to map one case onto another. Afterwards,
the two dataframes are compared on a cell to cell basis. If a cell in the 'OS_time' column matches, a marker for OS_time increases by one. The same is true for the 'OS' column. Once all 
cells have been compared, the 'OS_time_compare' is compared with 'OS_compare' and the length of 'survival'. If all three are equal, the program is successful, otherwise, the program fails. 
A list of missing cases, their submitter_id, OS_time, and OS is also printed. 
'''
def comparison(survival, xena, column):
    os_time_compare = 0 #used to track successful OS_time comparisons
    os_compare = 0 #used to track successful OS comparisons.
    missing_list = [] # cases that are missing from the xena file are appended to this list.
    for cell in range(len(survival.index)): # cell is a marker that goes through each row of 'survival'. 
        search = survival.loc[cell, column[2]] # search is a string that corresponds to specified row's submitter_id
        if (xena[column[2]].eq(search)).any() == True: # check if that submitter_id exists in the 'xena' dataframe. If so. we compare the other two columns. 
            row, col = np.where(xena == search) #finds row and column where submitter_id is located.
            xena_time = xena.loc[row[0], column[0]] #xena_time corresponds to the OS_time at specified row and column named 'OS_time'
            survival_time = survival.loc[cell, column[0]] # survival_time corresponds to the OS_time value at row 'cell'.
            xena_status = xena.loc[row[0], column[1]] # xena_status is the OS value at the specified row (based on submitter_id)
            survival_status = survival.loc[cell,column[1]] # survival_status corresponds to the OS value at row 'cell'.
            if xena_time == survival_time: # check if the two OS_times are equivalant. 
                os_time_compare += 1
            if xena_status == survival_status: # check if the two OSs are equivalent.
                os_compare += 1
            print(search) 
        else: # if 'search' (submitter_id) does not exist, then the case is missing from the xena file
            print("\nCase is missing from survival data matrix (submitter_id):")
            print(search)
            print("\n" + column[0] + " value:")
            print(survival.loc[cell, column[0]])
            print("\n" + column[1] + " value:")
            print(survival.loc[cell, column[1]])
            case_info = [search, survival.loc[cell, column[0]], survival.loc[cell, column[1]]] # the survival data on the case is then saved to 'missing_list'.
            missing_list.append(case_info)
    if os_compare == len(survival.index) and os_time_compare: # Checks if Os_time_compare and OS_compare are equivalent to the length of the 'survival' dataframe.
        print("\nsuccess!")
        print("number of cases compared: ")
        print(os_compare)
    else: # if they are not equivalent the program fails. 
        print("\nFailed")
        print("\nCases successfully compared:" + str(os_compare) + "/" + str(len(survival.index)))
        print("\nList of missing cases [Submitter_id, OS_time, OS]:")
        print(missing_list)
#gets the submitter_ids of a project using a post request.
#submitter_id_list is a list of submitter_ids. 
submitter_id_list = getSubmitterId(project, submitter_id_keys, submitter_id_fields, submitter_id_filters)

# gets all other survival data using another post request. 
# This includes vital_status (OS), days_to_death, days_to_last_follow_up, etc. (OS_time)
getData(submitter_id_list, survival_keys, survival_fields, survival_filter)

# all of this data is then formatted into a Dataframe
survival_df = formatData(column_names)
print(survival_df)
# the xena file TSV is also formatted into a dataframe.
xena_df = xenaFormat(xena_file)
print(xena_df)
# The two dataframes are compared and any missing cases or incorrect data is printed out at the end. 
comparison(survival_df, xena_df, column_names)
