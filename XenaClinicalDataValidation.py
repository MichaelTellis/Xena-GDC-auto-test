import requests
import json
import sys
import pandas as pd
import numpy as np
import sys

'''
	Xena clincal data validation is a script that compares clinical data that the Xena ETL code imported and clinical data 
imported independently from the GDC. The program first reads the Xena TSV file in order to get the relevant fields that
will be retrieved from the GDC (getFields). It then flips each field if it has a '.' in it, since the field names in the
Xena file are in opposite order (flipFields). The field: 'samples.submitter_id' is missing and needs to be added afterward
(addSubmitterId). Then, the xena file itself is converted to a dataframe (xenaFormat). The filter used by this script's 
request is 'samples.submitter_id', so we need 'samples.submitter_id' from the xena dataframe which will then be used as 
values to include in the request (getFilter). With all of the preperation complete, a request is sent to the cases endpoint
of the GDC, of all clinical data relevant to the fields recieved from 'getFields', and filtered by the 'samples.submitter_id'
from "getFilter". This request returns a json response which is then unpeeled to remove uneccessary dictionaries and lists
(unpeel). Ultimately, a list is returned by (getData). Then, the data is compared in "compareData()". Every element in 
the list response from getData() corresponds to a single case. So compareData() loops through each element. First, it 
format the data in "formatData()" to a dictionary, where each key corresponds to a field, and each value corresponding to 
the relevant data of that field. Some cases have multiple samples, so only those without 'ffpe scrolls' or 'blood derived normal',
are included in the data_dict. The formatted data_dict is passed to treatments(), a function where it checks the length of 
all values in data_dict with a key with 'treatements' in it, and then extends these values to the length of the longest list,
which is what happens in the xena file. The function submitterId() saves each 'samples.submitter_id' to a list. searches for
the specific 'samples.submitter_id' for the case we are in, by searching for the value of data_dict['samples.submitter_id'].
This submitter_id will then be used to find the spcific row in the xena dataframe with this samples.submitter_id. This tells
us the row with data relevant to this sample. Now with the row, we can go through each key in data_dict and match the value
with a cell in the xena dataframe. If the cell is a float value, it means that the value is a number or a pd Nan. And we compare
this float value with the value in the data_dict. If they match, then it passes and we move on to the next value. The other 
data type found in the xena dataframe is a string, which can either a normal string, or an empty string. We compare this
to the value found in data_dict. If they match, then it passes and we move on to the next cell. Once every value of every case is 
compared. We check the saples compared with the original samples found in the xena file. This is since some cases have multiple
samples, but we only end up comparing one. If there are missing samples, they are saved to a list ('missing') that will be sent through
the program a second time.  
'''
############################## GLOBAL #############################################
# list of all fields that will be requested from the GDC
fields = []  

# dataframe with all clinical data from the xena_file
xena_df = pd.DataFrame()

# list of 'samples.submitter_id' in the xena_file
filter_list = []

# List of all data retrieved from the GDC request
gdc_data = []
# Any samples that were not compared in the first run of program, and will be sent in a second time.
missing = []
# Any samples that were found in the GDC retrieved data, but not in the xena_file
missing_samples = []
############################### CONSTANT ##########################################

# Xena TSV file name. 
xena_file = "/Users/michaeltellis/Downloads/CTSP-DLBCL1/Xena_Matrices/CTSP-DLBCL1.clinical.tsv"

if len(sys.argv) != 2:
	print("incorrect arguments")
	sys.exit(0)
xena_file = sys.argv[1] #xena file path and name, second argument


#keywords used throughout program
keyword = ["treatments","samples.submitter_id",
 "submitter_id", "samples", "sample_type", 
 "Blood Derived Normal", "FFPE Scrolls", "submitter_id.samples"]


'''
unpeels unneccessary dictionaries and lists returned by GDC API. resJson is the GDC API response in Json format. However,
after being unpeeled it is a list. 
'''
def unpeel(resJson): # json format
    resJson = json.loads(resJson)
    resJson = resJson.get("data")
    resJson = resJson.get("hits")
    return resJson # list. 



'''
getFields reads the first line of 'xena_file' and splits it by \t, since xena_file is a TSV. then each element in 'line'
is stripped, and appended to 'fields', which now has the list of all fields in the xena_file.
'''
def getFields(xena_file):
	file = open(xena_file, "r+")
	line = file.readline() 
	file.close()
	line = line.split("\t")
	fieldList = []
	for i in line:
		fieldList.append(i.strip())
	return fieldList



'''
Due to the fields being flipped about a '.' in the xena_file, they need to be flipped again in order to return to the
format the GDC uses. This is done by taking in two arguments, the list of fields. and the cvhracter we are flipping about (arg).
'''
def flipFields(field_list, arg):
	flipped_field_list = []
	# for loop runs through each field in the field list. 
	for field in field_list:
		flipped = arg.join(field.split(arg)[::-1]) # splits them by arg, and flips them around, before joining them together again.
		flipped_field_list.append(flipped) # this new flipped field is appnded to 'flipped_field_list'.
	return flipped_field_list 



'''
addSubmitterId() adds 'samples.submitter_id' to the list of fields. 
'''
def addSubmitterId(fields):
	fields.append(keyword[1]) #keyword list is found in constants
	return fields



'''
xenaFormat turns the xena_file into a dataframe called xena_df. 
'''
def xenaFormat(file):
	xena_format = pd.read_csv(file, sep= '\t') 
	return xena_format


'''
getFilter() turns the column named 'submitter_id.samples' into a list of all submitter_ids found
in the Xena_file.
'''
def getFilter(xena_df):
	filters = xena_df[keyword[7]].tolist()
	return filters



'''
getData() takes in two arguments. The lsit of fields called 'fields' and the list of values for the filter, called 'filter_list'.
fields is the list of fields read from the xena_file, while filter list is a list of submitter_ids read from the xena dataframe.
A post request is sent to the GDC cases_endpt asking for a json response. This json response is then unpeeled into a list which 
is what the function returns. This list has all clinical data of the fields in 'fields' from the samples in 'filter_list'.
'''
def getData(fields, filter_list):
	fields = ",".join(fields) # the list of fields needs to be a comma seperated string when sent to the GDC

	cases_endpt = "https://api.gdc.cancer.gov/cases" # cases endpoint

	# The filter used is 'samples.submitter_id'. while 'value' is the submitter ids included in the request. 
	filters = {
		"op": "in",
		"content":{
			"field": "samples.submitter_id",
			"value": filter_list
			}
	}


	# A POST is used, so the filter parameters can be passed directly as a Dict object.
	params = {
    	"filters": filters,
		"fields": fields,
		"format": "json", #request for a json response
    	"size": "20000"
    	}

	# The parameters are passed to 'json' rather than 'params' in this case
	response = requests.post(cases_endpt, headers = {"Content-Type": "application/json"}, json = params)

	# turns the response into a readable json called 'responseJson'
	responseJson = json.dumps(response.json(),indent=2)
	#unpeels 'responseJson' into a list
	responseJson = unpeel(responseJson)
	return responseJson




'''
FormatData loops through the a dictionary with multiple nested dictionaries. This dictionary corresponds to data from a single case.
If the value in the key value pair is not a list or dictionary, it is added to a dicitonary called data_dict. The key is the field name
while the value is the data from the GDC, corresponding to the field. Sometimes, fields pop up multiple times in the dictionary and would
lead to duplicate keys. This is prevented by having the data from these identitical fields to be appended to the same key value pair, with
the data appended to the value which is a list. If the type of the value is a list or dictionary, we loop through each element
and check if it is not a list of dicitonary and add it to data_dict. Once every key, value pair is added to data_dict, it is returned. 
'''
def formatData(i, filter_list):
	data_dict = {}
	# each key is a field 	
	for key, value in i.items():
		# if the type is a string, float, None, or int it is added to data_dict with the key being the field name
		if type(value) != list and type(value) != dict:
			# if key already exists in data_dict, just append it to data_dict's value
			if key in data_dict:
				data_dict[key].append(value)
			# if not, create new key, value pair
			else:
				data_dict[key] = [value]
		# if the type is a list, we have to loop through every element in the list
		elif type(value) == list:
			for j in value:
				# If the type of j is not a list or dict, add it to data_dict
				if type(j) != list and type(j) != dict:
					if key in data_dict:
						data_dict[key].append(j)
					else: 
						data_dict[key] = [j] 
				# if the type of j is a dict, loop through the dict
				elif type(j) == dict:
					for name, item in j.items():
						# This checks if j corresponds to the 'samples' field group. As it has keyword[4] ('sample_type')
						# as one of the items. 
						if keyword[4] in j:
							#Now if 'sample_type' is in j, we check if the value of 'sample_type' is 'FFPE scrolls' or 'Blood derived normal'/
							#If they are, break the loop as they are not included in xena data
							if j[keyword[4]] == keyword[5] or j[keyword[4]] == keyword[6]:
								break
							# if the samples'submitter_id is not found in 'filter_list' which is the list of submitter_ids, then also break the loop.
							elif j[keyword[2]] not in filter_list:
								# However, this sample is saved to missing_samples, which is printed out at the end. 
								missing_samples.append(j[keyword[2]])
								break
						# if the type of the item is not a list or dict it is added to data_dict
						if type(item) != list and type(item) != dict:
							if key + '.' + name in data_dict:
								data_dict[key + '.' + name].append(item)
							else:
								data_dict[key + '.' + name] = [item]
						# if the type of the item is a list, we have to loop through one more time
						elif type(item) == list: 
							for k in item:
								# if the type of k is a dict, we loop through it as well. 
								if type(k) == dict:
									for marker, element in k.items():
										# if the type of the element is not a list or dict it is added to data_dict
										if type(element) != list and type(element) != dict:
											if key + '.' + name + '.' + marker in data_dict:
												data_dict[key + '.' + name + '.' + marker].append(element)
											else: 
												data_dict[key + '.' + name + '.' + marker] = [element]
										# If it is, an error will be printed out
										else:
											print("ERROR missing data path found in list, dict, list, dict, layer 5")
								# if the type of k is not a list or dict, we add it to data_dict
								elif type(k) != list and type(k) != dict:
									if key + '.' + name in data_dict:
										data_dict[key + '.' + name].append(item)
									else:
										data_dict[key + '.' + name] = [item]
								# If it is none of these, then an error is printed out
								else: 
									print("ERROR missing data path found in list, dict, list, layer 4")
						# if it is none of these an error is printed out
						else:
							print("ERROR missing data path found in list, dict, layer 3")
							print(name)
							print(item)
				# if it is none of these an error is printed out.
				else:
					print("ERROR missing data path found in list, layer 2")
		# if the type of the value is a dict, then we loop through each name, item pair in the dict
		elif type(value) == dict:
			for name, item in value.items():
				# if the type of the item is not a list or dict, we add it to data_dict
				if type(item) != list and type(item) != dict:
					if key + '.' + name in data_dict:
						data_dict[key + '.' + name].append(item)
					else:
						data_dict[key + '.' + name] = [item]
				# if  the type of the item is a list, we loop through the list
				elif type(item) == list:
					for j in item:
						# if the type of j is not a list or a dict, we add it to data_dict
						if type(j) != list and type(j) != dict:
							if key + "." + name in data_dict:
								data_dict[key + '.' + name].append(j)
							else:
								data_dict[key + '.' + name] = [j]      
						# if it is, then an error is printed out               
						else: 
							print("ERROR missing data path found in dict, list, layer 3")
				# if the type of the item is a dict then we loop through the marker, element pair of the dict.
				elif type(item) == dict:
					for marker, element in item.items():
						# if the type of the element is not a list or dict, it is added to data_dict
						if type(element) != list and type(element) != dict:
							if key + '.' + name + '.' + marker in data_dict:
								data_dict[key + '.' + name + '.' + marker].append(element)
							else:
								data_dict[key + '.' + name + '.' + marker] = [element]
						# if it is, an error is printed out
						else:
							print("ERROR missing data path found in dict, dict, layer 3")
				# if it is none of the above an error is printed out.
				else: 
					print("ERROR missing data path found in dict, layer 2")
		# if it is none of the above an error is printed out
		else:
			print("ERROR missing data path found, layer 1")
	# dictionary with all clinical data organized is returned. 
	return data_dict


'''
treatments() handles the treatments edge case. The xena file makes sure that the length of each field with 'treatments' is the number of treatments of the specific
case. So if there are 5 treatments, the number of elements in the list would be 5. In order to do that in data_dict, we find the maximum length of all values
where the key is 'treatments' and extend all list with 'treatments' to that length with empty strings. 
'''

def treatments(data_dict): # data_dict is dictionary with all clinical data.
	#treatment_list will include all keys with 'treatments'
	treatment_list = []
	for key ,value in data_dict.items(): 
		# check every key in data_dict for 'treatments', if 'treatments' is in the key, add it to treatment_list.
		if keyword[0] in key:
			treatment_list.append(value)
	# Finds length of each list, and saves it on list_len
	list_len = [len(i) for i in treatment_list]
	#If no treatments, then the length of list_len will be 0.
	if len(list_len) > 0:
		#list_len now int value.
		list_len = max(list_len)
		for key, value in data_dict.items():
			if keyword[0] in key:
				#If the length of the data_dict list is less then list_len, extend it by the difference. 
				if len(data_dict[key]) < list_len:
					extend = list_len - len(data_dict[key])
					data_dict[key].extend(["\'\'"] * extend)
	return data_dict


'''
submitterId() saves all samples.submitter_id(s) found in data_dict, and appends them to submitter_ids. Since this function is called
within a for loop, 'submitter_ids' wil have all samples.submitter_id(s) saved in a list. 
'''
def submitterId(data_dict, submitter_ids): # data_dict is dictionary with all clinical data. submitter_ids is list of all samples.submitter_ids. 

	for key, value in data_dict.items():
		# checks if key is equal to 'samples.submitter_id', if so append the str form of the value. 
		if keyword[1] == key:
			submitter_ids.append(str(value[0]))
			break
	return submitter_ids


'''
searchSubmitter_id() saves the 'submitter_id' to the search variable in string form. 'submitter_id' is different from 'samples.submitter_id'
since it is the name of the case, not sample. We will use the 'submitter_id' in compareXena() to find the row that correspondes to the data
in data_dict.
'''
def searchSubmitterId(data_dict): # data_dict is dictionary with all clinical data of a specific case.
	search = None
	for key, value in data_dict.items():
		# if the key is equal to 'submitter_id', then search is equal to the value.
		if keyword[2] == key:
			search = value
			search = str(search[0])
			break
	return search



'''
compareXena() compares the the values inside data_dict with the corresponding xena cell and checks if they match. First, the keys 
must be reversed. Then we find the value at the row and column of xena_df, which is assigned to 'xena'. Next, check the type of xena, 
which is either float, string or None. If Xena is float type, we compare it to each element in value (which is list type). If they're
equivlant, total_comp (total num of comparisons) and comp_count (number of equivlant comparisons) both increase by 1. Additionally, 
the value in xen is replaced by an empty string in case there are duplicate values. Otherwise, only total_comp increases by one. 
An important edge case is when xena is empty, but float type it is actually a pandas Nan, so we check if value is None or an empty
 string. The same is done if xena is string or None type. However, since all list values in the xena_file are converted to strings
 when converting to a dataframe, empty strings become '""' instead of '', if xena has '""', we check if value is None or '""'. Once
 every key, value pair is compared, check if total_comp is equivlant to comp_count. If so, the function returns True (successful 
 comparison) or False (unsuccesful).
'''
def compareXena(data_dict, search): # data_dict is dictionary with all clinical data of a specific case. Search is the submitter_id of this case.
	#total number of comparisons
	total_comp = 0 
	#total number of successful comparisons.
	comp_count = 0
	#checks the position when search is in the xena dataFrame.
	row, col = np.where(xena_df == search)
	#saves the row as an int value. 
	row = int(row[0])

	# loop through key, value pair in data_dict.
	for key, value in data_dict.items():
		# key is flipped about '.' in order to match the xena_df.
		key = ".".join(key.split(".")[::-1])
		#xena is the value at the specific row equal to the value of 'row' and column equal to the value of 'key'.
		xena = xena_df.loc[row,key]

		# If xena is float type:
		if type(xena) == float:
			# Loop through each element in value. 
			for i in range(len(value)):
				# if key has 'samples' then check if the Xena value is equal to the data_dict value. 
				if keyword[3] in key:
					if xena == value[0]:
						# if it is, remove the value. 
						xena = xena.replace(value[0], '', 1)
						# increase both markers by 1. 
						total_comp += 1
						comp_count += 1
				# if xena is equivlant to the element in value.
				elif xena == value[i]:
					# remove the value from xena
					xena = xena.replace(value[i], '', 1)
					# increase both markers by 1.
					total_comp += 1
					comp_count += 1
				# if xena is a pandas Nan, and this element in value is None type or '""'.
				elif pd.isna(xena) and (value[i] == None or value[i] == "\'\'"):
					# increase both markers by 1.
					total_comp += 1
					comp_count += 1
				# if xena is none of the above, only increase total comparisons. 
				#Print out the key, value pair and type, as well as xena. 
				else: 
					total_comp += 1
					print(key)
					print(value[i])
					print(type(value[i]))
					print(xena)
					print(type(xena))
					break
					break
					break
		# if xena is string type:
		elif type(xena) == str:
			# loop through each element in value
			for i in range(len(value)):
				# if key has 'samples' then check if the Xena value is equal to the data_dict value. 
				if keyword[3] in key:
					if xena == value[0]:
						# if it is, remove the value from xena.
						xena = xena.replace(value[0], '', 1)
						#both markers increase by 1.
						total_comp += 1
						comp_count += 1
				# if xena is an "empty string"('" "')(Since everything in lists are converted to a string, the quotations are also included)
				#and if value is none or '""'
				elif (value[i] == None or value[i] == "\'\'") and ("\' \'" in xena):
					# remove '" "' from xena. 
					xena = xena.replace("\' \'", '', 1)
					#both markers increase by 1
					total_comp += 1
					comp_count += 1
				#if the string form of value[i] is in xena.
				elif str(value[i]) in xena:
					#While in data_dict numbers are in int form, in the xena file, they are floats.
					if type(value[i]) == int:
						#change value to a float. 
						value[i] = float(value[i])
					# Then remove the str(value) from xena
					xena = xena.replace(str(value[i]), '', 1)
					#both martkers increase by 1. 
					total_comp += 1
					comp_count += 1
				# if xena is none of the above, only increase total comparisons. 
				#Print out the key, value pair and type, as well as xena. 
				else: 
					total_comp += 1
					print(key)
					print(value[i])
					print(type(value[i]))
					print(xena)
					print(type(xena))
					break
					break
					break
		# if xena is None type:
		elif type(xena) == None:
			# loop through each element in value
			for i in range(len(value)):
				# if the value is equal to xena:
				if value[i] == xena:
					#both markers increase by 1
					total_comp += 1
					comp_count += 1
				# if xena is none of the above, only increase total comparisons. 
				#Print out the key value pair and type, as well as xena.
				else:
					total_comp += 1
					print(key)
					print(value)
					print(type(value[i]))
					print(xena)
					print(type(xena))
					break
					break
					break
	# if the total number of comparisons is equal to the number of succesful comparisons. return True (successful)
	if total_comp == comp_count:
		return True
	# if not, return False (unsuccessful)
	else: 
		return False



'''
compareData() calls the functions formatData(), treatments(), submitterId(), searchSubmitterId(), compareXena() in a loop in order
to run through every single case in responseJson. It creates a data_dict, a dictionary with all clinical data in formatData(). Runs
data_dict through treatments() in order to handle the treatments edge case (expained in treatments()). Saves all 'samples.submitter_id's 
in submitter_id by running submitterId(). Gets the submitter_id of the current case by running searchSubmitterId(). Then calls
compare Xena to compare the xena_df with data_dict. If the comparison is successful, sample_compare (number of successful sample
comparisons) increases by one. If the program fails, the sample is appended to failed_samples. num_sample_compare (total number of samples
compared) increases by 1. At the end, if any samples were not compared, they are saved to missing_data and then run through the whole
program a second time. 
'''
def compareData(responseJson, filter_list): # responseJson (probably a bad nameas it is a list) is a list with all clinical data from the selected samples.
	#filter_list is a list of all 'samples.submitter_ids' in the xena file. 

	#submitter_ids will have all of the 'samples.submitter_ids' that were compared.
	submitter_ids = []

	#num_sample is the length of filter_list and is just the number of samples. 
	num_sample = len(filter_list)

	#num_sample_comapre is number of samples compared.
	num_sample_compare = 0
	#sample_compare is number of samples successfully compared
	sample_compare = 0
	#failed_cases is the list of 'submitter_id' where the comparison was unsuccessful.
	failed_cases = []
	#missing_data  is the list of samples that were not compared.
	missing_data = []

	#loops througheach element in responseJson, i corresponds to a single case.
	for i in responseJson:
		#data_dict is a dictionary of all data in i. 
		data_dict = formatData(i, filter_list)

		#handles the treatments edge case.
		data_dict = treatments(data_dict)

		#gets the 'samples.submitter_id' from data_dict
		submitter_ids = submitterId(data_dict, submitter_ids)

		#gets the 'submitter_id' from data_dict. Search will be used in compareXena
		search = searchSubmitterId(data_dict)

		#compares data_dict with xena_df. 
		#If the function returns true, it means the comparison was successful. 
		if compareXena(data_dict, search) == True:
			#number of successfully compared samples increases by 1.
			sample_compare +=1
		#otherwise, the comparison was unsucessful.
		else:
			#the 'submitter_id' of the case is saved to failed_cases
			print("failed")
			failed_cases.append(search)
		#total number of samples compared increases by 1.
		num_sample_compare +=1			


	print(sample_compare)
	print(num_sample_compare)
	print(num_sample)
	print(submitter_ids)
	sub = set(submitter_ids)
	# checks if x (samples.submitter_id) is not in the list of subitter_ids compared. If it isn't it is added to missing_data.
	missing_data = [x for x in filter_list if x not in sub]
	print("Any missing sample data that was not compared: ")
	print(missing_data)
	print("Any cases that failed the comparison:")
	print(failed_cases)
	return missing_data


#fields is a list of fields in the xena_file
fields = getFields(xena_file)

#flips the fields around about '.'
fields = flipFields(fields, '.')

#adds 'samples.submitter_id' to the list of fields
fields = addSubmitterId(fields)

#converts xena_file to dataFrame called xena_df
xena_df = xenaFormat(xena_file)

#list of "samples.submitter_id" that will be used as filters when retrieving data from the GDC.
filter_list = getFilter(xena_df)

#Retrieves relevant clinical data from the GDC
gdc_data = getData(fields, filter_list)

#compares the GDC data to the xena dataFrame. 
#any samples that were missing from the comparison are added to 'missing'.
missing = compareData(gdc_data, filter_list)

#checks if there are any missing samples.
if len(missing) > 0:
	#If so retrieve data a second time and compare a second time.
	gdc_data = getData(fields, missing)
	missing = compareData(gdc_data, missing)
# Any samples that are missing and not in the list of submitter_ids from the xena file, must be missing from the xena file. 
for sample in missing_samples:
	if sample not in filter_list:
		print(sample)

