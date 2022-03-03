#!/usr/bin/env python3
"""
Author: Aura Zelco

Title: buscoParser.py

- General description:
This script retrieves the sequences of orthologs found with BUSCO in all target 
species, and creates a FASTA-like file for each BUSCO hit, with one sequence 
per species (top ortholog in the BUSCO output). 

- Procedure:
1. checks if the directories exist
2. imports the BUSCO ID list: these are IDs that are present in all species
3. imports the species name in a list, which have to match the name of the BUSCO subfolders
4. creates a nested dictionary: for each species (keys of outer dictionary), an 
inner dictionary is created containing the gene ids (keys) and sequences (values) 
from the .faa files
5. retrieves sequences of orthologs for each BUSCO id and writes a new output file 

- Usage:
This script retrieves the sequences for unique BUSCO orthologs and outputs one 
file for each BUSCO hit, containing the orthologoues sequences, one per species.

It is run in the command line as: python buscoParser.py -busco_id BUSCO_ID 
-species SPECIES_ID -busco_path BUSCO_PATH -faa_path PROTEIN_SEQ_PATH 
[-output_path OUTPUT_PATH]

-List of user-defined functions:
1. dir_path: checks if directory exists
2. retrieve_ids: from a txt file, cretaes a list with every line as an element
3. species_dictionaries: creates a nested dictionary containing the species, 
the gene ids and the sequences from the .faa files
4. find_files: checks if file exists
5. find_top_ortho: finds the top ortholog gene id from BUSCO outputs
6. retrieve_sequences: retrieves the sequence of the ortholog from the species_dictionaries

- List of imported modules:
1. argparse: a module which is used to input the different parameters. 
2. os: module to check if given directory exists and to retrieve filenames
3. re: regex search used here to calculate the GC content

- Possible bugs:
1. The script only checks if the directories and files exist. There are no 
further checks on the files content in the directories. 

"""
#%% IMPORT MODULES

import argparse
import os
from os import path
import re

#%% ARGPARSE

# description of program, printed when -h is called in the command line
usage = 'This script retrieves the sequences for unique BUSCO orthologs and outputs one file for each BUSCO hit, containing the orthologoues sequences, one per species. '

# creates an ArgumentParser object which has assigned name 'parser'
parser = argparse.ArgumentParser(description=usage)


# specification of input file
parser.add_argument(
    '-busco_id',
    metavar = 'BUSCO_ID',
    dest = 'busco_id',
    type=argparse.FileType('r'), # readable file
    required=True,   # needs to be inserted in the command line
    help="BUSCO unique IDs" 
    )

parser.add_argument(
    '-species',
    metavar = 'SPECIES_ID',
    dest = 'species',
    required=True,  # needs to be inserted in the command line
    type=argparse.FileType('r'), # readable file
    help="input species id, have to match the name of the BUSCO folders" 
    )

parser.add_argument(
    '-busco_path', 
    metavar = 'BUSCO_PATH',
    dest = 'busco_path',
    type=os.path.abspath, # extracts the absolute path, easier to navigate through the tree
    required=True,  # needs to be inserted in the command line
    help="path to input directory" 
    )

parser.add_argument(
    '-faa_path', 
    metavar = 'PROTEIN_SEQ_PATH',
    dest = 'faa_path',
    type=os.path.abspath, # extracts the absolute path, easier to navigate through the tree
    required=True,  # needs to be inserted in the command line
    help="path to protein sequences directory" 
    )

parser.add_argument(
    '-output_path', 
    metavar = 'OUTPUT_PATH',
    dest = 'out_path',
    type=os.path.abspath, # extracts the absolute path, easier to navigate through the tree
    default = os.path.curdir, # the default is the present working directory
    help="optional path to output directory" 
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
    # if the directory exists, return the direcotry path as it is
    if os.path.isdir(string):
        return string
    # if not, raise an error
    else:
        raise argparse.ArgumentTypeError("{path_dir} is not a valid path")

'''RETRIEVE_IDS

    Parameters
    ----------
    filename : string
        file path as string

    Returns
    -------
    id_list: list
        list containing all lines present in the file

'''

# retrieves the sequences for a given file
def retrieve_ids(filename):
    # initiates a list
    id_list=[]
    # reads the files
    for line in filename:
        # appends the line to the list
        id_list.append(line.strip())
    # returns list
    return(id_list)

'''SPECIES_DICTIONARIES

    Parameters
    ----------
    species_list : list
        list containing the species names
    faa_dir: path
        folder absolute path where .faa are found

    Returns
    -------
    species_faa_dict: nested dictionary
        outer dictionary contains species as keys,
        values are inner dictionaries with faa gene ids (keys) and sequences (values)

'''

# creates a nested dictionary
def species_dictionaries (species_list, faa_dir):
    # initializes the outer dictionary
    species_faa_dict={}
    # loops thorugh the species
    for sp in species_list:
        # creates the inner dictionary
        sp_dict={}
        # retrieves the faa file path
        faa_filepath = '{}/{}.faa' .format(faa_dir, sp)
        # opens the file
        with open(faa_filepath) as faa_file:
            # reads each line
            for line in faa_file:
                # if we have a header
                if line.startswith('>'):
                    # extract the gene name (first element in tab-delimited file), while removing the >
                    sp_key = line.split('\t')[0].strip('>')
                    # retrieves the sequence in the next line and removes the newline
                    sp_dict[sp_key] = next(faa_file).strip()
        # adds the inner dictionary to the outer dictionary
        species_faa_dict[sp] = sp_dict
    # returns the nested dictionary
    return species_faa_dict
    
'''FIND_FILES

    Parameters
    ----------
    filename : string
        file path as string

    Returns
    -------
    busco_id_path: string
        file absolute path

'''

# checks if BUSCO files exists
def find_files(busco_id, busco_folder, sp):
    # retrieves the full path to the BUSCO output files
    busco_id_path=f'{busco_folder}/{sp}/run_apicomplexa_odb10/hmmer_output/initial_run_results/{busco_id}.out'             
    # if the file exists
    if path.exists(busco_id_path):
        # returns the absolute path
        return(busco_id_path)
    # if not
    else:
        # raises an error
        raise ('The file {} does not exist.' .format(busco_id_path))
    
'''FIND_TOP_ORTHO

    Parameters
    ----------
    filename : string
        file path as string

    Returns
    -------
    genename: string
        gene id of top ortholog

''' 

# finds the top ortholog gene id
def find_top_ortho(filename):
    # set flag to False
    first_hit = False
    # opens the file
    with open(filename) as gene_file:
        # reads each line
        for line in gene_file:
            # if first_hit is False
            if not first_hit:
                # if the line contains _g, which is always present in the gene name
                if bool(re.search('_g', line))==True:
                    # saves the genename
                    genename = line.split()[0]
                    # sets the flag to True so genename is saved only for the first gene
                    first_hit=True
    # returns the gene id
    return(genename)

            
'''RETRIEVE_SEQUENCES

    Parameters
    ----------
    id_list: list
        BUSCO unique id list
    busco_folder: string
        absolute path to BUSCO result folder
    species_list: list
        species list
    species_dict: nested dictionary
        outer dictionary contains species as keys,
        values are inner dictionaries with faa gene ids (keys) and sequences (values)

    Returns
    -------
    fileOut: output file
        output file containing a FASTA-like format of species, ortholog id and sequence
        one per BUSCO id

'''

# retrieves the sequences and outputs the FASTA-like files
def retrieve_sequences (id_list, busco_folder, species_list,species_dict):
    # loops through the BUSCO id list
    for element in id_list:
        # opens a new output file
        fileOut = open((f'{out_path}/{element}.parser.output'), 'w')
        # loops through the species
        for sp in species_list:
            # retrieves the filename from the BUSCO output folder
            sp_filename=find_files(element, busco_folder, sp)
            # retrieves the gene/protein name
            protname_sp = find_top_ortho(sp_filename)
            # retrieves the sequence of the ortholog
            prot_seq = species_dict[sp][protname_sp]
            # write the species, the gene name and the sequence to the output file
            fileOut.write('>{}\tgene_name={}\n{}\n' .format(sp, protname_sp, prot_seq))
        # closes the output
        fileOut.close()
                

#%% MAIN

# checks if BUSCO directory exists
busco_dir=dir_path(args.busco_path)
# checks if faa files directory exists
faa_dir=dir_path(args.faa_path)
# checks if output directory exists
out_path=dir_path(args.out_path)


# creates the BUSCO id list
id_list = retrieve_ids(args.busco_id)
# stores the species name
species_list=retrieve_ids(args.species)
# creates one dictionary with nested species dictionaries
species_dict=species_dictionaries(species_list, faa_dir)


# retrieves the sequences for each BUSCO id, and ouputs a file for each
retrieve_sequences(id_list, busco_dir, species_list, species_dict)

