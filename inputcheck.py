""" inputcheck module of TreeDomPy.

This module contains the funtions to check query and gbk files.
"""

#=========
# Modules
#=========

from Bio import SeqIO

def is_fasta(query):
    """Check if query is in fasta format
    Args: query
    Returns: True or False
    """
    with open(query, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        isfasta=any(fasta) # False when `fasta` is empty, i.e. wasn't a FASTA file

        return isfasta  

def is_gbk(gbk_list):
    """Check if gbks are in gbk format
    Args: list of gbks
    Returns: True or False
    """
    isgbk_dict=dict()
    for gbk in gbk_list:
    	with open(gbk, "r") as handle:
    		gbk_TF = SeqIO.parse(handle, "genbank")
    		isgbk_dict[gbk]=any(gbk_TF) # False when `gbk_TK` is empty, i.e. wasn't a GBK file
    return isgbk_dict  

