#!/usr/bin/env python3
"""
Author: Aura Zelco

Title: removeBird.py

- General description:
This script removes the scaffolds obtained from the datParser.py from the clean genome obtained from removeScaffold.py. 

- Procedure:
1. creates a set from the scaffold file from the command line
2. reads the genome file (infile) and extracts the ids and sequences, line by line
3. writes to the output file only the scaffolds and sequences not contained in the scaffolds set

- Usage:
This program removes certain scaffolds based on the id and outputs a new filtered genome file

It is run in the command line as: python removeBird.py -i input_genome [-o output] [-s scaffolds.txt]

- List of imported modules:
1. argparse: a module which is used to input the different parameters. 

- Possible bugs:
1. If no output file is given, the results will be written in the default output file,
whether that filename exists already in the directory. Thus, there is a risk of overwriting. 

"""
#%% IMPORT MODULES
import argparse

#%% ARGPARSE

# description of program, printed when -h is called in the command line
usage = 'This program removes certain scaffolds based on the id and outputs a new filtered genome file.'

# creates an ArgumentParser object which has assigned name 'parser'
parser = argparse.ArgumentParser(description=usage)

# specification of input file
parser.add_argument(
    '-i',
    metavar = 'INFILE',
    dest = 'infile',
    type=argparse.FileType('r'),# readable file
    required=True,  # needs to be inserted in the command line
    help="input genome" 
    )

# adds arg which writes an output: here the output is optional, so if not present the output will go to default name
parser.add_argument(
    '-o',
    dest = 'outfile',
    metavar = 'OUTFILE ',
    type=argparse.FileType('w'), # writable file
    default='filtered_genome.output',
    help='specify output file name (default is filtered_genome.output)'
    )

# adds arg which is required, containg the IDs of the scaffolds to be removed
parser.add_argument(
    '-s', 
    dest = 'scaffolds',
    metavar = 'SCAFFOLDS',
    default='bird_scaffolds.txt',
    type=argparse.FileType('r'), # readable file
    help="scaffolds to be removed"
    )


# returns result of parsing 'parser' to the class args
args = parser.parse_args()

#%% MAIN

# creates a list from the scaffold file
scaffoldlist = []
for line in args.scaffolds:
    scaffoldlist.append('>'+line.strip())

# reads the genome input file
for line in args.infile:
    # selects the headers
    if line.startswith('>'):
        # removes the newline
        scaffold_id = line.strip()
        # checks if the scaffolds is one of those to be removed
        if scaffold_id not in scaffoldlist:
            # if not, stores the sequence (found in the next line) in seq and removes the newline
            seq = next(args.infile).strip()
            # prints the ids and sequences to the output file
            print('{}\n{}' .format(scaffold_id, seq), file=args.outfile)

