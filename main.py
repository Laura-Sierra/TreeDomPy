# ===========================================================================================================================
# title           :main
# description     :TreeDomPy is an user-friendly tool to make blastP alignments, phylogenetic trees and to annotate domains.
# development     :In python using a Tkinter graphical interface.
# author          :Laura Sierra Heras
# last_revised    :16/01/2020
# version         :v.0.1
# usage           :python main
# ===========================================================================================================================


#=========
# Modules
#=========

from os.path import normpath, basename
import os

import tkinter as tk
import tkinter.messagebox
from tkinter import ttk
from tkinter import filedialog
from PIL import Image, ImageTk
from Bio import SeqIO
import argparse

import inputcheck
import gbkparser
import blaster
import muscle
import prosite
import plotter
import organizer

#==========
# Argparse
#==========
# TreeDomPy do not need arguments, but argparse is used to do --help | -h
parser=argparse.ArgumentParser(description="Description: TreeDomPy is an user-friendly tool \
to (1) align queries vs the proteins of genebank(s) using blastP, \
(2) make trees after muscle alignement \
and (3) detect domains using Posite database. \
It is developed in Python using a Tkinter graphical interface. \
For more information about what kind of files you need and about output files. \
Please, use TreeDomPy (python main.py) and press on help button.")
parser.parse_args()

#=================
# Tkinter Window
#=================
def main():
    class Root(tk.Tk):
        def __init__(self):
            super(Root, self).__init__()

            # Window style definition: tittle, size and background
            self.title("TreeDomPy")
            self.minsize(400, 400)
            self.img = ImageTk.PhotoImage(Image.open('window_img.png'))
            self.panel = ttk.Label(self, image=self.img)
            self.panel.grid(columnspan= 2, row=8)
            self.configure(bg='SeaGreen1')
     
            # Button to choose gbk files. Settings:
            self.gbk_labelFrame = ttk.LabelFrame(self, text = "Choose genbank  Files")
            self.gbk_labelFrame.grid(column = 0, row = 1, padx = 20, pady = 20)
            self.gbk_button = ttk.Button(self.gbk_labelFrame, text = "Browse...",command = self.gbk_chooser)
            self.gbk_button.grid(column = 1, row = 1)
     
            # Button to choose query file. Settings:
            self.query_labelFrame = ttk.LabelFrame(self, text = "Choose query.fa File")
            self.query_labelFrame.grid(column = 1, row = 1, padx = 20, pady = 20)
            self.query_button = ttk.Button(self.query_labelFrame, text = "Browse...",command = self.query_chooser)
            self.query_button.grid(column = 1, row = 1)

            # Slider and button to choose coverage cut-off. Settings:
            self.cov_labelFrame = ttk.LabelFrame(self, text = "Coverage")
            self.cov_labelFrame.grid(column = 0, row = 4, padx = 20, pady = 20)
            self.cov_var= tk.IntVar()
            self.cov_scale = tk.Scale(self, from_=0, to=100,variable=self.cov_var)
            self.cov_scale.grid(column = 0, row = 3)
            tk.Button(self.cov_labelFrame, text="Click to Save", command=self.check_cov_scale).grid(column = 0, row = 4)

            # Slider and button to choose identity cut-off. Settings:
            self.id_labelFrame = ttk.LabelFrame(self, text = "Identity")
            self.id_labelFrame.grid(column = 1, row = 4, padx = 20, pady = 20)
            self.id_var = tk.IntVar()
            self.id_scale = tk.Scale(self, from_=0, to=100,variable=self.id_var)
            self.id_scale.grid(column = 1, row = 3)
            tk.Button(self.id_labelFrame, text="Click to Save", command=self.check_id_scale).grid(column = 1, row = 4)

            # Help button. Settings:
            self.help = tkinter.Button(self, text = "Help", command = self.help)
            self.help.grid(column = 0, row = 5)

            # Run button. Settings:
            self.run = tkinter.Button(self, text = "Run", command = self.run)
            self.run.grid(column = 1, row = 5)  


        # Procedure to choose gbk files and to show file names selected
        def gbk_chooser(self):
            self.gbk_filename = filedialog.askopenfilenames(parent=self,title='Choose the genbank files')
            gbk_list=root.tk.splitlist(self.gbk_filename)
            self.gbk_list=gbk_list
            gbk_txt=""
            for gbk in gbk_list:
            	gbk=basename(normpath(gbk))
            	gbk_txt=gbk_txt+str(gbk)+"\n"
            self.gbk_label = ttk.Label(self.gbk_labelFrame, text = "")
            self.gbk_label.grid(column = 1, row = 2)
            self.gbk_label.configure(text = gbk_txt)

        # Procedure to choose query file and to show file name selected
        def query_chooser(self):
            self.query_filename = filedialog.askopenfilename(parent=self,title='Choose the file with your query/queries')
            query_list=root.tk.splitlist(self.query_filename)
            self.query_list=query_list
            try:
                query_txt=str(basename(normpath(query_list[0])))+"\n"
                self.query_label = ttk.Label(self.query_labelFrame, text = "")
                self.query_label.grid(column = 1, row = 2)
                self.query_label.configure(text = query_txt)
            except:
                pass

        # Procedure to save coverage cut-off chosen
        def check_cov_scale(self):
            self.cov_label = ttk.Label(self.cov_labelFrame, text = "")
            self.cov_label.grid(column = 0, row = 2)
            self.cov_var_get=self.cov_var.get()
            self.cov_label.configure(text = self.cov_var_get)
    	
        # Procedure to save identity cut-off chosen
        def check_id_scale(self):
            self.id_label = ttk.Label(self.id_labelFrame, text = "")
            self.id_label.grid(column = 1, row = 2)
            self.id_var_get=self.id_var.get()
            self.id_label.configure(text = self.id_var_get)	

        #==============
        # HEP MESSAGE:
        #==============
        # Pep8 sais that when we split lines they should be indent as the first one, but in this case it is imposible.
        # Help button takes indentation as real spaces
        def help(self):
            help_text="""TreeDomPy NEEDS:
1. A file or a set of genbankfiles. To choose several files use Ctrl+Alt.
2. A fasta file with one or several queries.
3. A coverage cut-off. You must press button save, once you choose it with the slider.
4. An identity cut-off. You must press button save, once you choose it  with the slider.
Files used are copied to Data+Date directory

TreeDomPy GENERATES:
1. gbk_multi.fasta: fasta file from your gbk(s).
2. blast_hits.tsv: blast with your preferences.
3. blast_hits.fasta: fasta of the hits.
4. domains_hits.txt: hits domains and its explanation.
5. muscle_results.fasta: hits alignment with muscle.
6. tree.nwx: NJ-tree.
Files generated are placed in Results+Date directory.

TreeDomPy FIGURES:
Blast, tree and domains."""
            tkinter.messagebox.showinfo("Help", help_text)

        # Run
        def run(self):
        #=============
        # CHECK INPUT:
        #=============
            # Check prosite.dat user has unzip prosite.dat.zip
            if not os.path.exists("prosite.dat"):
                warning_massage="Please, unzip prosite.dat.zip. \nLeave prosite.dat file in the same folder where main.py can be found"
                tkinter.messagebox.showwarning("Warning!",warning_massage)
                return
                
            # Check gbks, query, id and cov    
            warning_dict={'gbk': ["Please, choose at least one gbk file.\n\n","The following file(s) is/are not in gbk format:\n"], 
                        'query': ["Please, choose the fasta file with your query/queries.\n\n","The file with your query is not in fasta format.\n\n"], 
                        'cov': 'Please, choose a coverage cut-off.\nRemember to save it once you choose it.\n\n',
                        'id': 'Please, choose an indentity cut-off.\nRemember to save it once you choose it.\n\n'}
            warning_massage=""           
            # Check gbk(s)
            try:
                if self.gbk_filename == tuple():
                    raise Exception
                isgbk_dict=inputcheck.is_gbk(self.gbk_filename)
                wrong_files=""
                for key, value in isgbk_dict.items():
                    if not value:
                        wrong_files=wrong_files+"\t· "+basename(normpath(key))+"\n"               
                if wrong_files!="":
                    warning_massage=warning_massage+warning_dict.get("gbk")[1]+wrong_files+"\n"
            except:
                warning_massage=warning_massage+warning_dict.get("gbk")[0]
            # Check query
            try:
                isfasta=inputcheck.is_fasta(self.query_filename)
                if not isfasta:
                    warning_massage=warning_massage+warning_dict.get("query")[1]
            except:
                warning_massage=warning_massage+warning_dict.get("query")[0]
            # Check coverage
            try:
                existe=self.cov_var_get
            except:
                warning_massage=warning_massage+warning_dict.get("cov")
            # Check identity
            try:
                existe=self.id_var_get
            except:
                warning_massage=warning_massage+warning_dict.get("id")
            if warning_massage!="":
                tkinter.messagebox.showwarning("Warning!",warning_massage+ "If you need more information.\nPlease, click on help button.")
                return 
            else:
                tkinter.messagebox.showinfo("Let´s go!","Your inputs have been validated.\nPress OK to start the process.\nIt may take a while, PLEASE BE PATIENT.\n\nRemember:\n-Save and close figures to allow the porgram continue.\n-The more files and queries you use, the longer it will take.")                




        #==========================
        # RUN IF INPUT IS CORRECT 
        #==========================  
            try:
                gbk_multi_fasta=gbkparser.parse_gbk(self.gbk_filename) # Pass from gbk format to fasta
                cont=1 # Loop to go through queries
                for seq_record in SeqIO.parse(self.query_filename, "fasta"):
                    query= "query.fasta" # This file will contain one query at a time (the next in each iteration)
                    output_handle = open(query, "w")
                    if str(seq_record.id) !="":
                        output_handle.write(">"+str(seq_record.id)+"\n") # print identifier of the hit
                    else:
                        output_handle.write(">Query\n") # print >Query in case header is empty
                    output_handle.write(str(seq_record.seq)+"\n") # print sequence of the hit
                    output_handle.close()
                    hits_sseqid_list,dict_plot_blast=blaster.blast(self.cov_var_get,self.id_var_get) # Make blast
                    if type(hits_sseqid_list) is str:
                        tkinter.messagebox.showerror("Error", "Opps! Something went wrong.\n"+hits_sseqid_list) # Show blast error to user in case the is one
                        os.remove("query.fasta") # We remove files created so that they cannot gum up the future work of the user
                        os.remove("gbk_multi.fasta")
                        return
                    plotter.make_blast(dict_plot_blast) # Plot Blast
                    hits_fasta=blaster.hits_to_fasta("query.fasta",gbk_multi_fasta,hits_sseqid_list) # Pass results of blast to fasta in order to do the tree and look for domains
                    dict_to_plot,max_seq_len,pattern_list=prosite.domain_parser() # Search domains
                    muscle.muscle(hits_fasta) # Make tree
                    if os.path.exists("tree.nwx"):
                        plotter.make_tree(self.gbk_list) # Plot tree
                    else: # Muscle only generates tree file if there are at least 3 sequences, in case they are less we tell the user why we cannot show the tree.
                        tkinter.messagebox.showwarning("Warning", "Zero hits detected.\nIt is not possible to make a tree of one sequence.\nPress OK to see the domains of the sequences.")
                    plotter.make_domains(dict_to_plot,max_seq_len,pattern_list) # Plot Domains
                    
                    # Organize data and results files 
                    if cont==1: #For query1 its a little bit diferent. We create results and data folders
                        time_var=organizer.results_folder(cont,0)
                    else:
                        organizer.results_folder(cont,time_var)
                    cont=cont+1
                os.remove("query.fasta") # We remove that provisional file. That now contains the last query.
                organizer.data_folder(self.query_list,self.gbk_list,time_var)
            except: #Presumably this should not happend
                tkinter.messagebox.showerror("Error", "Opps! Something went wrong.\nRevise your files and make sure you have installed blast and muscle.\nWe do not know exactly why the error ocurred.\nSo please, delete any intermediate file created in the package folder so that they cannot gum up your future work.")
    root = Root()
    root.mainloop()


if __name__ =='__main__':
    main()

