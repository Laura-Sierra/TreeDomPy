""" muscle module of TreeDomPy.

This module is used for muscle alignment.
"""

#=========
# Modules
#=========

from subprocess import Popen, PIPE,call


def muscle(hits_fasta):
	"""Muscle alignement and tree.nwx generation 
    Args: hits in fasta format
    """

    # Alignement
	muscle_process = Popen(['muscle','-in',hits_fasta, "-out", "muscle_result.fasta"], stdout=PIPE, stderr=PIPE)
	muscle_error = muscle_process.stderr.read().decode('utf-8')
	muscle_result = muscle_process.stdout.read().decode('utf-8')
	muscle_process.stderr.close()
	muscle_process.stdout.close()

	# Tree
	muscle_process = Popen(['muscle','-maketree','-in',"muscle_result.fasta", "-out", "tree.nwx", "-cluster", "neighborjoining"], stdout=PIPE, stderr=PIPE)
	muscle_error = muscle_process.stderr.read().decode('utf-8')
	muscle_result = muscle_process.stdout.read().decode('utf-8')
	muscle_process.stderr.close()
	muscle_process.stdout.close()


	return

	    
