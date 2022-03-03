#!/usr/bin/env python3
"""
Author: Aura Zelco

Title: removeScaffold.py

- General description:
This script takes an input a genome file, and filters out scaffolds according 
to their length and GC content. The output is saved in a new file. 

- Procedure:
1. defines a function to calculate the GC content in a sequence
2. opens input genome
3. if line is header, stores it in seq_id
4. if line is a sequence, concatenates the lines in variable seq
5. when a new id is found, checks if the sequence has the necessary requirements (GC and minimum length)
6. if so, writes the id and the sequence to the output file

- Usage:
This program takes an input genome file and filters the scaffolds according 
to the GC content and minimum length parameters

It is run in the command line as: python removeScaffold.py -i input_genome [-o output] [-gc GC percentage] [-min_length minimum scaffoldlength] 


- List of functions:
1. GC_cont: calculates the GC percentage in a sequence

- List of imported modules:
1. argparse: a module which is used to input the different parameters. 
2. re: a module which allows for use of regular expressions (hereafter RegEx)

- Possible bugs:
1. This script does not allow for a lower threshold of GC content; thus, the 
scaffolds can be filtered only if GC is lower than a value.
2. The GC content is calculated on the whole length of the sequence, including 
of N are present; therefore, the calculated GC content might be lower than the 
actual one, and more sequences are filtered.  
3. If no output file is given, the results will be written in the default output file,
whether that filename exists already in the directory. Thus, there is a risk of overwriting. 
"""
#%% IMPORT MODULES
import argparse
import re

#%% ARGPARSE 
# description of program, printed when -h is called in the command line
usage = 'This program takes an input genome file and filters the scaffolds according to the GC content and minimum length parameters.'

# creates an ArgumentParser object which has assigned name 'parser'
parser = argparse.ArgumentParser(description=usage)

# specification of input genome file
parser.add_argument(
    '-i',
    metavar = 'INFILE',
    dest = 'infile',
    type=argparse.FileType('r'), # readable file
    required=True, # needs to be inserted in the command line
    help="input genome" 
    )

# adds arg which writes an output: here the output is optional, so if not present the output will go to default name
parser.add_argument(
    '-o',
    dest = 'outfile',
    metavar = 'OUTFILE ',
    type=argparse.FileType('w'), # writable file
    default='genome.output',
    help='specify output file name (default is genome.output)'
    )

# arg which will store the optional GC content threshold
parser.add_argument(
    '-gc_cont', 
    dest = 'gc_thresh',
    metavar = 'GC_THRESH',
    type=int,
    default=30,
    help="optional GC threshold in percentage (default 30)"
    )

# arg to store the minimum length needed for the sequences
parser.add_argument(
    '-min_length', 
    dest = 'min_length',
    metavar = 'MIN_LENGTH',
    type=int,
    default=3000,
    help="optional minimum scaffold length (default 3000)"
    )

# returns result of parsing 'parser' to the class args
args = parser.parse_args()

#%% USER-DEFINED FUNCTIONS
''' GC_cont

    Parameters
    ----------
    seq : string
        nucleotide sequence

    Returns
    -------
    gc_cont: float
        GC percentage in the sequence
'''

# defines a new function to count the GC content
def GC_cont(seq):
    # counts all matches of G or C in the string
    gc = len(re.findall('[GC]', seq))
    # calculates the percentage using the whole sequence length
    perc_gc = gc*100/len(seq)
    # returns the percentage
    return(perc_gc)
#%% MAIN
# intializes an empty string
seq = ''

# reads the lines in the input genome file
for line in args.infile:    
    # selects the header lines    
    if line.startswith('>'):
        # takes just the scaffold name
        seq_id = line.strip('\n').split()[0]
        # if there is a sequence
        if seq:
            # if the sequence is longer than 3000 bases
            if len(seq) >= args.min_length:
                # calculates the GC
                gc_seq = GC_cont(seq)
                # if the GC is lower than the given threshold
                if gc_seq <= args.gc_thresh: 
                    # print the ID, GC content and sequence to the file
                    print('{}\tGC={:.2f}\n{}' .format(seq, gc_seq, seq_id), file=args.outfile)
                    # re-initializes the sequence to empty
                    seq = ''
        # for the first ID, there is no sequence
        else:
            # so we just print the first header
            print('{}' .format(seq_id), file=args.outfile)
    # if it is a sequence line
    else:
        # removes the newline and makes all characters upper case
        line = line.strip('\n').upper()
        # adds to the string, since in this case the sequence is on multiple lines
        seq += line
# when we reach the end of the file, we still have one sequence to write to the output file
if seq:
    print('{}' .format(seq), file=args.outfile)
            


