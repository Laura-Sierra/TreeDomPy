""" prosite module of TreeDomPy.

This module is intended to search for domains. It used prosite database.

I know that the best option would have been to make a file with the patterns just like python searchs.
So that we do not have to make pattern modification each time.
But in this way, I can show how the conversion is made (prosite pattern -> re pattern).
And I do not modify prosite database, in case the user wants to add something to the program.
"""

#=========
# Modules
#=========

import re

from Bio.ExPASy import Prosite,Prodoc
from Bio import SeqIO
from prettytable import PrettyTable

def domain_parser():
	"""Search domains of each hit
    Returns: dictionary with values useful to graphic, len of the hit with the largest sequence, list of detected patterns
    """
	input_handle="blast_hits.fasta"
	# File with: (1) the domains of each hit and (2) domains information
	domains= "domains_hits.txt"
	output_handle = open(domains, "w")
	output_handle.write("#This file contains the domains of each hit.\n#At the bottom, you will find detail information of all the domains detected.\n")
	output_handle.write("#We strongly recommend to open this file with Visual Studio Code.\n#Because when the names of the domains are too large, in regular editors the table looks awful.\n")
	output_handle.write("#Here it is only showed how many times a pattern is present.\n#In the figure of the domains you will find the position of each domain.\n\n")
	accession_list=[] 	# List of prosite.doc accessions of the domains that had been found
	domains_dict=dict() # dictionary that saves matches
	count=1
	max_seq_len=0 # Keep larger sequence to plot x-axe
	# Loop to go through hits
	for seq_record in SeqIO.parse(input_handle, "fasta"):
		output_handle.write(str(seq_record.id)+"\n") # print identifier of the hit
		output_handle.write(str(seq_record.seq)+"\n") # print sequence of the hit
		if len(seq_record.seq)>max_seq_len:
			max_seq_len=len(seq_record.seq)
		# Make a table for each hit with the domains, that contains the following fields: name, accession, description and pattern
		x=PrettyTable()
		x.field_names=["name","accession","description","pattern","repetitions"]

		# Loop to go through prosite domains
		handle = open("prosite.dat","r")
		records = Prosite.parse(handle)
		for record in records:
			# prosite.dat preparation for parsing
			# {} -> [^]
			pattern = record.pattern.upper()
			pattern = pattern.replace("{", "[^")
			pattern = pattern.replace("}", "]")	
			# - -> ""
			pattern = pattern.replace("-", "")	
			# . -> ""
			pattern = pattern.replace(".", "")	
			# X|x -> "[ARNDCQEGHILKMFPSTWYV]"
			AAS="[ARNDCQEGHILKMFPSTWYV]"
			pattern = pattern.replace("x", AAS)
			pattern = pattern.replace("X", AAS)	
			# () -> {}
			pattern = pattern.replace("(", "{")
			pattern = pattern.replace(")", "}")	

			# >] -> ]?$
			pattern = pattern.replace(">]", "]?$")	

			#  <  -> ^
			#  >  -> $
			pattern = pattern.replace("<", "^")	
			pattern = pattern.replace(">", "$")	
			if pattern != "":
				# Look if the hit contains the current patter
				if re.search(r""+str(pattern), str(seq_record.seq).upper()): # if found
					if record.pdoc not in accession_list:
						# Save pdoc accession in the list of prosite.doc accessions
						# if it is not already
						accession_list.append(record.pdoc)
					matches = re.finditer(r""+str(pattern), str(seq_record.seq).upper())
					reps=0
					for match in matches: # save all matches in a dictionary to plot them later
						domains_dict[count]=[seq_record.id, len(seq_record.seq),record.name,match.start(),match.end()]
						count=count+1
						reps=reps+1
					x.add_row([record.name,record.accession,record.description,record.pattern, reps]) # add found domain to table

		output_handle.write(str(x)+"\n") # add table of hit to domains_hits.txt

	# At the end of the tables, print information of all the domains that had been found
	output_handle.write("\n")
	record_text_list=DocParser(accession_list)
	for text in record_text_list:
		output_handle.write(text)
	return (domains_dict,max_seq_len,accession_list)


def DocParser(accession_list):
	"""Get record.txt info for domains detected
    Args: list of accesions of the domains found
    Returns: record.txt of all patterns
    """
	handle = open("prosite.txt")
	records = Prodoc.parse(handle)
	record_text_list=[]
	try:
		# Loop to go through prosite.doc entries. 
		for record in records:
			if record.accession in accession_list: # If an entry is in the list of the already found domains
				record_text_list.append(record.text) # Save it
	except:
		print(False)
	return record_text_list


