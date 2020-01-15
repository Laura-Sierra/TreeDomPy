""" organizer module of TreeDomPy.

This module is to organaize input and output files 
"""

#=========
# Modules
#=========

import os
import shutil
import datetime
from os.path import normpath, basename


def results_folder(cont,time_var):
	"""Organize output data
    Args: numer of query we are testing (as the first one works different), time to add to the folder name
    """
	# try/except to check if a folder with this name already exists
	try:
		# Create Results and Data folders.
		# Date added to folder, to avoid problem of existence
		if cont==1: # When we are with the secord query, we do not need to create Data and results
			now = datetime.datetime.now()
			time_var=now.strftime("%Y-%m-%d_%H:%M")
			os.mkdir("Results_"+time_var,0o777)
			os.mkdir("Data_"+time_var,0o777)
		os.mkdir("Results_"+time_var+"/query"+str(cont),0o777)


	except: 
		# Presumably, this will never happend, but just in case.
		print("The folders called Data_"+time_var+" and Results_"+time_var+" already exists.")

	# Moving files to Results
	my_current_directory = os.getcwd()
	shutil.move(my_current_directory+"/blast_hits.tsv", my_current_directory+"/Results_"+time_var+"/query"+str(cont)+"/blast_hits.tsv")
	try:
		shutil.move(my_current_directory+"/tree.nwx", my_current_directory+"/Results_"+time_var+"/query"+str(cont)+"/tree.nwx")
	except:
		pass
	shutil.move(my_current_directory+"/blast_hits.fasta", my_current_directory+"/Results_"+time_var+"/query"+str(cont)+"/blast_hits.fasta")
	shutil.move(my_current_directory+"/domains_hits.txt", my_current_directory+"/Results_"+time_var+"/query"+str(cont)+"/domains_hits.txt")
	shutil.move(my_current_directory+"/muscle_result.fasta", my_current_directory+"/Results_"+time_var+"/query"+str(cont)+"/muscle_result.fasta")

	return time_var


def data_folder(query_list,gbk_list,time_var):
	"""Organize input data and save gbk_multi.fasta
    Args: query file, gbk files, time to add to the folder name
    """
	my_current_directory = os.getcwd()
	# Moving files to Data
	for gbk in gbk_list:
		gbk_base=basename(normpath(gbk))
		shutil.copyfile(str(gbk), my_current_directory+"/Data_"+time_var+"/"+gbk_base)
	for n_query in range(len(query_list)):
		query_base=basename(normpath(query_list[0]))
		shutil.copyfile(str(query_list[0]), my_current_directory+"/Data_"+time_var+"/"+query_base)
	# Moving gbk_multi.fasta to Results as it is required for all the queries
	shutil.move(my_current_directory+"/gbk_multi.fasta", my_current_directory+"/Results_"+time_var+"/gbk_multi.fasta")
	return


