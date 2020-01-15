""" gbkparser module of TreeDomPy.

This module is used to extract proteins from gbks
"""

#=========
# Modules
#=========

from Bio import SeqIO

def parse_gbk(gbk_files):	
	"""Pass gbk to fasta (extract proteins)
    Args: gbk files
    Returns: multifasta of the gbks
    """
	gbk_multi_fasta = "gbk_multi.fasta"
	output_handle = open(gbk_multi_fasta, "w")
	# Loop to go through gbk files
	gbk_list=list(gbk_files)
	for file in gbk_list:
		input_handle  = open(file, "r")
		for seq_record in SeqIO.parse(input_handle, "genbank"):
			# Loop to go through all the proteins of gbk file
		    for seq_feature in seq_record.features :
		        if seq_feature.type=="CDS":
		        	try:
		        		# Fasta file format:
		        		# > protein_id/of/organism_id product [source]
		        		# sequence
				        output_handle.write(">%s %s#%s\n%s\n" % (
				        seq_feature.qualifiers['protein_id'][0],
				        #seq_record.id,
				        seq_feature.qualifiers['product'][0].replace(" ","-").replace(";","-").replace("(","").replace(")","").replace(",",""),
				        seq_record.annotations["source"].replace(" ","-").replace(";","-").replace("(","").replace(")","").replace(",",""),
				        seq_feature.qualifiers['translation'][0]))
		        	except:
		        		pass

		input_handle.close()

	output_handle.close()
	return gbk_multi_fasta


