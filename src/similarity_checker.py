from Bio import SeqIO
from Bio import Align
import os

import scipy as sc
import pandas as pd

import Bio.Phylo
# To self: need to tidy up code into functions, possibly a 'Mystery' Class for the unknown sequence that allows comparison, probability table, and phylogeny functions? 

# read in the dog breed database
database_file = ".\data\dog_breeds.fa"
database_list = list(SeqIO.parse(database_file, "fasta"))

dog_database = {}
for dog in database_list: # when optimising: Bio.SeqIO.to_dict(), and Bio.index() might get 1-step dictionary?
    description_list = dog.description.split('[')
    breed = description_list[7][6:-2]
    dog_database[breed, dog.id] = dog.seq

# read in the mystery sequence
sequence_directory = ".\data\seq"
for filename in os.listdir(sequence_directory): # generic to allow any named fasta file (tests x2: is there a file in this folder? is this a fasta file?)
    if filename.endswith(".fa"):
        path = os.path.join(sequence_directory, filename)
        mystery_seq = SeqIO.read(path, "fasta")

# comparison module
# tests: length of sequence = the same? 
alignment_scoring = {}
for k, v in dog_database.items():
    print("Now checking alignment with:", k[0])
    aligner = Align.PairwiseAligner()
    alignment = aligner.align(v, mystery_seq.seq)
    alignment_scoring[k] = alignment.score
print(max(dict, key = dict.get)) # to make this output nicer eg an information file in the results folder