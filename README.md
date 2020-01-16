# TreeDomPy
TreeDomPy is an user-friendly tool to make blastP alignments, phylogenetic trees and to annotate domains using Python and a Tkinter graphical interface.

![alt text](https://raw.githubusercontent.com/Laura-Sierra/TreeDomPy/images/0_window.png)
Note: If you feel as though the program is not working. Please, be patient, some processes take a while to compute. In case there is a problem, an error window will immediately appear.
## Requirements
Python version: Python 3.x.

Modules: Tkinter, numpy, Bio, argparse, matplotlib, pylab and prettytable.

External programs: Blast and Muscle.

Prosite DataBase: As I had problems with prosite files provided I uploaded mine. You must unzip prosite.dat.zip. It was too big to upload it as a normal file.

prosite.doc is changed to prosite.txt and I delete some kind of \r at the end of a line, just that.

I know that the best way to proceed would have been to make a file with the patterns as required by the re module (since there are slight changes with those of prosite database) so as not to make those changes every time I use the domain parser. But I though that in this way, I can show how the conversion is made (prosite pattern -> re pattern).
And I provide the original prosite database, in case the user wants to add something to the program.

## Usage

```sh
python main.py
```
TreeDomPy can be run in UNIX and MAC systems with python 3.x. GUI is optimized for Tkinter in UNIX, weird aspects for main Tkinter window may be observed in MacOSX systems.

## Input files
TreeDomPy requieres genbank or several genbank files. It also requieres a fasta or multifasta file containing all the queries. The input files are loaded using the Browse button. Remember format of input files will be checked.

This files will be saved in the Data_TimeOfExperiment Folder.

This program has been validated with several gbks and queries of:
  - Bacteria (Escherichia coli, Salmonella enterica, Shigella dysenteriae and Yersinia pestis). This four gbks are provided with a query of example (TestData Folder)
  - Animals (Mus musculus)
  - Plants (Arabidopsis thaliana).

## Parameters
TreeDomPy only needs two parameters: coverage and idenity cut-off values for blastP (e-value is set to 0.000001). They are chosen with a slider bar, once chosen they must be saved. 

## Output
TreeDomPy generates several files stored in the Results_TimeOfExperiment folder (in case there are several queries, they are ordered by query)
1. gbk_multi.fasta: fasta file from your gbk(s).
2. blast_hits.tsv: blast results with your preferences of identity and coverage.
3. blast_hits.fasta: fasta of the hits.
4. domains_hits.txt: hits domains and its explanation.
5. muscle_results.fasta: hits alignment with muscle.
6. tree.nwx: NJ-tree.

TreeDomPy also generates 3 figures.
Remember to save the figures you like. You must close figures so that the programa could follow its workflow.
1. BlastP Results Figure
  - The coverage is represented with bars and identity with a color scale (query is the red one).
  - Remember if there are so many hits that you cannot see a clear figure, choose higher values for coverage and identity.
![alt text](https://raw.githubusercontent.com/Laura-Sierra/TreeDomPy/images/1_blast.png)
2. Tree Figure

  - If you want the figure to publish it in a paper but you do not like the colors of the figure. Run the program again. The colors are randomly generated. Therefore, you will have different designs for the same tree-figure.
![alt text](https://raw.githubusercontent.com/Laura-Sierra/TreeDomPy/images/2_tree.png)

3. Domains Figure

  - If you want the figure to publish it in a paper but you do not like the colors of the figure. Run the program again. The colors are randomly generated. Therefore, you will have different designs for the same domains-figure.
![alt text](https://raw.githubusercontent.com/Laura-Sierra/TreeDomPy/images/3_domains.png)


## Author
LAURA SIERRA HERAS

For any questions, feel free to contact me.
