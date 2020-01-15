""" blaster module of TreeDomPy.

This module is used to make blast and to get results of blast in a useful form
"""

#=========
# Modules
#=========

from subprocess import Popen, PIPE,call
import re


def blast(cov_cutoff,id_cutoff):
    """Make blastP
    Args: coverage cut-off, identity cut-off
    Returns: hits identifiers, dictionary with values useful to graphic blast results
    """
    blast_process = Popen(['blastp','-query',"query.fasta",'-subject',"gbk_multi.fasta",'-qcov_hsp_perc',str(cov_cutoff),'-evalue', '0.000001', '-outfmt',"6 qseqid qcovs pident evalue sseqid sseq qstart qend" ], stdout=PIPE, stderr=PIPE)
    #e-value set to 0.000001
    blast_error = blast_process.stderr.read().decode('utf-8')
    hits_eval_cov = blast_process.stdout.read().decode('utf-8')
    blast_process.stderr.close()
    blast_process.stdout.close()

	# Filter by identity and save hits in a file.tsv
	# Return a list of hit identifiers
    if not blast_error:
        filtered_blast="blast_hits.tsv"
        output_handle  = open(filtered_blast,"w")
        output_handle.write("#Query_id\tCoverage\tIdentity\tE-value\tSubject_id\tSubject_seq\tAlign_stat\tAlign_end\n")
        hits_eval_cov_list=hits_eval_cov.split("\n")
        hits_sseqid_list=[]
        dict_plot_blast=dict()
        cont=0
        for hit in hits_eval_cov_list:
            hit_list=hit.split("\t")
            try:
                if float(hit_list[2]) >= id_cutoff:
                    output_handle.write(hit+"\n")
                    dict_plot_blast[cont]=[hit_list[4],hit_list[2],hit_list[6],hit_list[7]]
                    cont=cont+1
                    if hit_list[4] not in hits_sseqid_list:
                        hits_sseqid_list.append(hit_list[4])
            except:
                pass
		
        output_handle.close()
        return hits_sseqid_list,dict_plot_blast
    else: 
        return (blast_error,0)

def hits_to_fasta(query,gbk_multi_fasta,hits_sseqid_list):
	"""Pass hits to fasta format
    Args: query file, gbk now in fasta, list of hit identifiers
    Returns: file of hits in fasta format
    """
	hits_fasta = "blast_hits.fasta"
	output_handle = open(hits_fasta, "w")
	# Add query to fasta for muscle alignement
	with open(query) as f_query:
	    for line in f_query:
	        output_handle.write(line.replace("#",""))
	# Add hits to fasta for muscle alignement
	with open(gbk_multi_fasta) as f_hits:
	    header = f_hits.readline()
	    sequence=f_hits.readline()
	    while header:
		    for hit_sseqid in hits_sseqid_list:
		    	if re.search(r""+hit_sseqid, header):
		    		output_handle.write(header+sequence)
		    header = f_hits.readline()
		    sequence=f_hits.readline()
	output_handle.close()
	return hits_fasta
