#!/usr/bin/env python3
"""
Author: Aura Zelco

Title: GCcontent.py

- General description:
This script calculates the overall GC content in all genome files present in a directory. 

- Procedure:
1. checks if the directory exists
2. if so, for each filename retrieves the absolute path
3. reads each file, and stores the genome in one string
4. calculates the GC content on this string 
5. prints the result to the console

- Usage:
This script calculates the overall GC content in all genome files present in a directory.

It is run in the command line as: python GCcontent.py -path PATH_DIR

-List of user-defined functions:
1. dir_path: checks if directory exists
2. GC_cont: calculates the GC percentage in a sequence
3. process_file: formats genome in one line, and applies GC_cont
4. process_directory: runs process_file on all files in the directory

- List of imported modules:
1. argparse: a module which is used to input the different parameters. 
2. os: module to check if given directory exists and to retrieve filenames
3. re: regex search used here to calculate the GC content

- Possible bugs:
1. The script only checks if the directory exists. There are no further checks on the genome files in the directory. 


"""
#%% IMPORT MODULES
import argparse
import os
import re

#%% ARGPARSE

# description of program, printed when -h is called in the command line
usage = 'This program calculates the GC content for every genome in the given directory.'

# creates an ArgumentParser object which has assigned name 'parser'
parser = argparse.ArgumentParser(description=usage)

parser.add_argument(
    '-path', 
    metavar = 'PATH_DIR',
    dest = 'path_dir',
    required=True,  # needs to be inserted in the command line
    type=os.path.abspath,  # extracts the absolute path, easier to navigate through the tree
    help="path to input directory" 
    )


# returns result of parsing 'parser' to the class args
args = parser.parse_args()

#%% USER-DEFINED FUNCTIONS

'''DIR_PATH

    Parameters
    ----------
    seq : string
        folder path as string

    Returns
    -------
    string: folder path as string

'''

# checks if directory is valid
def dir_path(string):
    # if the directory exists, return the directory path as it is
    if os.path.isdir(string):
        return string
    # if not, raise an error
    else:
        raise argparse.ArgumentTypeError("{path_dir} is not a valid path")

''' GC_CONT

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
    # calculates the length of the sequence, excluding Ns
    tot = len(re.findall('[GCAT]', seq))
    # calculates the percentage using the whole sequence length
    perc_gc = gc*100/tot
    # returns the percentage
    return(perc_gc)


'''PROCESS_FILE

    Parameters
    ----------
    filename : string
        a file in the folder given from the command line

    Returns
    -------
    prints the filename and the GC content
'''

# opens the genome file and calculates the GC content, which is then printed to the console
def process_file(filename):
    # initializes an empty string where the whole genome will be stored
    seq = ''
    # opens the file
    with open(filename, 'r') as file_input:
        # reads the files
        for line in file_input:
            # excluding the headers
            if not line.startswith('>'):
                # adds sequence to seq variable without the newline
                seq += line.strip()
    # after the file has been read, it caculates the GC content on the whole genome
    gc_file = GC_cont(seq)
    # prints the filename and the GC content to the console
    print('{}\n{:.2f}%' .format(filename, gc_file))
    

'''PROCESS_DIRECTORY

    Parameters
    ----------
    path_dir : string
        folder given from the command line

    Returns
    -------
    runs functions for each file
    
'''

# takes the directory from the command line
def process_directory(path_dir):
    # for each file in the directory
    for filename in os.listdir(path_dir):
        # joins the directory path and the filename, so we obtain the absolute path of the file to run the GC function on
        process_file(os.path.join(path_dir, filename))

#%% MAIN

# checks if directory exists
path_dir=dir_path(args.path_dir)
# if so, runs the function
process_directory(path_dir)


