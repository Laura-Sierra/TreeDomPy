""" plotter module of TreeDomPy.

This module contains the funtions to plot blastp, tree and domains. 
"""

#=========
# Modules
#=========

from random import randint

from Bio import Phylo
from Bio import SeqIO
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors
import pylab
import numpy as np



def make_blast(dict_to_plot):
    """Plot blastP results, showing coverge of hits and identity of aligned parts based on a color-scale.
    Args: dictionary with the necessary data to plot
    """
    # Set figura size
    fig = plt.figure(figsize=(30, 30), dpi=100)

    # Set legend with a colormap from 0-100, which will corresponf to id %
    n=100
    cmap = plt.get_cmap("plasma", 100)
    norm= matplotlib.colors.BoundaryNorm(np.arange(0,n+1)-0.5, n)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb=plt.colorbar(sm, ticks=np.arange(0,101,10))
    cb.set_label('Identity (%)')

    # Show query bar
    input_handle="query.fasta"
    for seq_record in SeqIO.parse(input_handle, "fasta"):
        max_seq_len=len(str(seq_record.seq)) 
        l=plt.hlines(y=str(seq_record.id), xmin=0, xmax=max_seq_len,linewidth=6, color="r")

    # Draw aligned part of each hit, the color depends on its identity
    cont=0
    for key, value in dict_to_plot.items():
        cont=cont+1
        color=cmap(float(value[1])/100)
        if cont<10:
            sep=" # "
        else:
            sep=" #"
        l = plt.hlines(y=value[0]+sep+str(cont), xmin=int(value[2]), xmax=int(value[3]),linewidth=6, color=color)
        

    plt.axis([0, max_seq_len, -1, cont+1]) # axis: max-x = length of query max-y = number of proteins to plot   
    plt.title('BlastP Results',fontsize = 20) # Title
    #Draw Blast alignement
    plt.show()
    return

def make_tree(gbk_list):
    """Plot tree, showing source of the hit with different colors
    Args: the list of gbk (to see the different sources)
    """

	# Read tree file
    tree = Phylo.read("tree.nwx", "newick")
    org_dict=dict()     
    # Keys = organisms of gbks (source) & Values = hits of this organism
    count=0 # add number to each hit, in orden to distinguish hits with the same name.
    for  clade in tree.find_clades():
    	if clade.name == None:
    		clade.name=""
    	else:
        	count=count+1
        	clade_name_list=str(clade.name).split("#") # [0] -> hit and [1] -> source organisms
        	clade.name=clade_name_list[0].replace("-"," ")
        	clade.name=str(count)+" "+clade.name
        	clade.name=clade.name[:37] # The label tree is up to 37 characters
        	try: # add new entry to ord_dict
        		if clade_name_list[1] not in org_dict.keys(): # if organism not already in dict
        			org_dict[clade_name_list[1]]=[clade.name] # add hit in a new list
        		else: # if it is already
        			org_dict[clade_name_list[1]].append(clade.name) # add hit
        	except: #for the queries
        		try: # for the fist query
        			org_dict["Query"].append(clade.name)
        		except: # for the second/thrid/ ... query in there are several
        			org_dict["Query"]=[clade.name]

    colors = ["r"] # Red for query
    for i in range(len(gbk_list)): # Make as many colors as gbk files
        colors.append('#%06X' % randint(0, 0xFFFFFF))

    col_dict=dict()
    # Tip label color dictionary -> Keys = hit label in tree & Values = Color of label
    leg_dict=dict()
    # Legend dictionary -> Keys = organisms of gbks (source) & Values =  its color
    cont=0
    for key, value in org_dict.items():
        if key == "Query": # query goes red
            for val in value:
                col_dict[val]=colors[0]
                leg_dict[key.replace("-"," ")]=colors[0]
        else:
            cont=cont+1
            for clade in tree.find_clades():
                if clade.name in value:
                    col_dict[clade.name[:37]]=colors[cont]
                    leg_dict[key.replace("-"," ")]=colors[cont] 


	# set the size and font of the figure
    fig = plt.figure(figsize=(30, 30), dpi=100)
    axes = fig.add_subplot(1, 1, 1)
    matplotlib.rc('font', size=6,weight='bold')
    # set legend. with leg_dict
    matplotlib.rcParams["legend.loc"] ='upper left'
    markers = [plt.Line2D([0,0],[0,0],color=color, marker='o', linestyle='') for color in leg_dict.values()]
    plt.legend(markers, leg_dict.keys(), numpoints=1)
    # set title
    plt.title('Tree', fontsize = 20)

    # Draw tree
    Phylo.draw(tree, axes=axes, label_colors=col_dict)
    return


def make_domains(dict_to_plot,max_seq_len,pattern_list):
    """Plot domains, each with a different color
    Args: dictionary with the necessary data to plot, len of the maximum hit, list of patterns found
    """
    fig = plt.figure(figsize=(30, 30), dpi=100)
    colors = [] # Make as many colors as possible pattern names
    for i in range(len(pattern_list)+10):
        colors.append('#%06X' % randint(0, 0xFFFFFF))
    hit_list=[] # list of proteins identifiers
    domains_dict=dict()     # Legend dictionary -> Keys = domain names & Values =  its color
    count=-1
    # Loop to go through all the pattens that had been found (of all proteins) and asing a color
    for key, value in dict_to_plot.items():
        if value[2] not in domains_dict.keys():
            count=count+1
            domains_dict[value[2]]=colors[count] # if domain name not in legend dictionary add it with its color 
   
    # Loop to go through all the pattens that had been found (of all proteins) and draw them
    for key, value in dict_to_plot.items():
        if value[0] not in hit_list:
            hit_list.append(value[0]) # if protein not in hit list add it
            l = plt.hlines(y=value[0], xmin=0, xmax=value[1],linewidth=2, color="c") # draw protein sequence based on length

        # Loop to search which color goes with the patten examined
        for key_col, value_col in domains_dict.items():
            if key_col == value[2]:
                colour=value_col
        # and draw pattern
        l = plt.hlines(y=value[0], xmin=value[3], xmax=value[4],linewidth=6, color=colour)

    plt.axis([0, max_seq_len+10, -1, len(hit_list)]) # axis: max-x = length of largest sequence max-y = number of proteins to plot
    # set legend. with domains_dict
    markers = [plt.Line2D([0,0],[0,0],color=color, marker='o', linestyle='') for color in domains_dict.values()] 
    plt.legend(markers, domains_dict.keys(), numpoints=1,bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',ncol=4, mode="expand", borderaxespad=0.)
    # set title
    plt.xlabel("Domains",fontsize=20)
    # Draw Domains
    plt.show()
    return
