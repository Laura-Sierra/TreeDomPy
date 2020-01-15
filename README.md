# TreeDomPy
TreeDomPy is an user-friendly tool to make blastP alignments, phylogenetic trees and to annotate domains using a Tkinter graphical interface.

![alt text](https://raw.githubusercontent.com/Laura-Sierra/TreeDomPy/images/0_window.png)

## Requires
Python version: Python 3.x.

Modules: Tkinter, numpy, Bio, argparse, matplotlib, pylab and prettytable.

External programs: Blast and Muscle.

Prosite DataBase: As I had problems with prosite files provided I uploaded mine. 

prosite.doc is change to prosite.txt and I delete some kind of \r at the end of a line, just that.

I know that the best way to proceed would have been to make a file with the patterns as required by the re module (since there are slight changes to those of prosite database) and I am making those changes every time I use the domain parser, but I though that in this way, I can show how the conversion is made (prosite pattern -> re pattern).
And I do not modify prosite database, in case the user wants to add something to the program.

## Usage

```sh
python main.py
```
TreeDomPy can be run in UNIX and MAC systems with python 3.x. GUI is optimized for Tkinter in UNIX, weird aspects for main Tkinter window may be observed in MacOSX systems.

## Input files
TreeDomPy requieres genbank or several genbank files. It also requieres a fasta or multifasta file containing all the queries. The input files are loaded using the Browse button. Remember format of input files will be checked.

## Parameters
TreeDomPy only needs two parameters: coverage and idenity cut-off values for blastP (e-value is set to 0.000001). They are chosen with a slider bar, once chosen they must be save. 

## Output
