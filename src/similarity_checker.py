from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align
import math
import os
import pandas as pd
import numpy as np
import scipy.stats as sc
import random
from Bio.Phylo import TreeConstruction

# To self: need to tidy up code into functions, possibly a 'Mystery' Class for the unknown sequence that allows comparison, probability table, and phylogeny functions? 

# read in the dog breed database
database_file = ".\data\dog_breeds.fa"
database_list = list(SeqIO.parse(database_file, "fasta")) #Note: parses seqeucne attributes as a string

breeds = [] # initialise lists for input into pandas dataframe
sample_id = []
sequence = []

for dog in database_list: # when optimising: Bio.SeqIO.to_dict(), and Bio.index() might get 1-step dictionary?
    description_list = dog.description.split('[')
    breed = description_list[7][6:-2]
    breeds.append(breed) # adds data for each dog to the dataframe row by row
    sample_id.append(dog.id)
    sequence.append("".join(dog.seq))

dog_database = pd.DataFrame (
    {
        "breed":breeds,
        "sample_id":sample_id,
        "sequence":sequence,
        "align_score":0
    }
)

# read in the mystery sequence
sequence_directory = ".\data\seq"
for filename in os.listdir(sequence_directory): # generic to allow any named fasta file (tests x2: is there a file in this folder? is this a fasta file?)
    if filename.endswith(".fa"):
        path = os.path.join(sequence_directory, filename)
        mystery_seq = SeqIO.read(path, "fasta")

# comparison module
# set up the aligner
# tests: length of sequence = the same? 
alignment_scoring = []
aligner = Align.PairwiseAligner()
aligner.mode = "local" # set to global for the global alignment function. Local alignment is slower than global.

# set up Karlin-Altschul parameters


# perform alignment for the mystery sequence against all sequences in the database
for i in dog_database.index:
    print("Now checking alignment with:", i, dog_database.loc[i, "breed"]) # progress bar
    alignment = aligner.align(dog_database.loc[i, "sequence"], mystery_seq.seq) # Smith-Waterman local alignment
    dog_database.loc[i, "align_score"] = alignment.score # writes alignment score to dataframe, saves
    
    # For local alignments only: Karlin-Altschul Algorithm (E- and P-value for statistical quality of alignment, same algorithm as BLASTn)
    score = alignment.score 
    e_value = 0.041 * len(mystery_seq) * len(dog_database.index) * math.exp(-0.267 * score) # probability of getting a score more than the one you generate above
    p_value = 1 - math.exp(-e_value) # estimates p-value from e-value. Assumes a Gumbel extreme value distribution of sequences.
    print(e_value, p_value)
    dog_database.loc[i, "E_value"] = e_value # writes alignment score to dataframe, saves
    dog_database.loc[i, "p_value"] = p_value # writes alignment score to dataframe, saves

# output table of aligned scores sorted by alignment score
dog_database = dog_database.sort_values("align_score")
dog_database.to_csv("results/final/similarity_alignment.csv")

# probability scores
# scipy.stats import pearsonr (standard correlation coefficient)
# def pearsonr_pval(x,y):return pearsonr(x,y)[1]
# corr = df.corr(method=pearsonr_pval) - generates a correlation matrix for all sequences vs all sequences in the database

# phylogeny tree
# global multiple sequence alignment for database of sequences
# calculator = DistanceCalculator('identity') - identity is the scoring model, but can check attr of calculator to see other models available
# above.get_distance(aln) - gives a DistanceMatrix object
# constructor = DistanceTreeConstructor(calculator, 'nj') - uses nearest neighbour joining method
# tree = constructor.build_tree(aln)
# bootstrapping to get the best tree: 
# trees = bootstrap_trees(aln, 100, constructor, majority_consensus) - generates just the consensus tree
# draw() to get a matplotlib rooted tree image output