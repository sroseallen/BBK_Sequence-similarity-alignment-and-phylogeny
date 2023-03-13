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
database_file = ".\data\dog_breeds_short.fa"
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
# tests: length of sequence = the same? 
alignment_scoring = []

# set up the aligner
aligner = Align.PairwiseAligner()
aligner.mode = "local" # set to global for the global alignment function. Local alignment is slower than global.

# set up Karlin-Altschul parameters
# counts of nucleotides in the whole database
nt_count = {"A":0, "T":0, "G":0, "C":0}
for seq in dog_database["sequence"]:
    for nt in nt_count:
        nt_count[nt] += seq.count(nt)

# counts the number of transitions between each nucleotide in each sequence
trans_counts = {"A":{"A":0, "T":0, "G":0, "C":0},
                "T":{"A":0, "T":0, "G":0, "C":0},
                "G":{"A":0, "T":0, "G":0, "C":0},
                "C":{"A":0, "T":0, "G":0, "C":0}
                } # nested dictionary: allows each nucleotide to have a count of nuclotides that have occured after it in sequence
for seq in dog_database["sequence"]: 
    for i in range(len(seq)-1): # required to prevent key error
        nt1 = seq[i]
        nt2 = seq[i+1]
        if nt1 not in nt_count or nt2 not in nt_count:
            continue # accounts for any gaps in database (algorithm requires ungapped database sequences for E-value estimation)
        trans_counts[nt1][nt2] += 1

# generate probabilities for all bases and base transtions across dog_database
base_prob = {}
for nt1 in nt_count:
    base_prob[nt1] = {}
    for nt2 in nt_count:
        if nt_count[nt1] == 0:
            base_prob[nt1][nt2] = 0
        else:
            base_prob[nt1][nt2] = trans_counts[nt1][nt2] / nt_count[nt1]
M = sum([nt_count[nt] for nt in nt_count])
L = sum([nt_count[nt] * nt_count[nt2] * base_prob[nt][nt2] for nt in nt_count for nt2 in nt_count])

l = -1 / (M * (L - 1)) # lambda* value in Karlin-Altshuls' Theorem 1 - sensitivity of alignment to matches
k_value = l * M * M
length_seq = len(mystery_seq)

# shuffle the mystery sequence several times and store (for p-value calculation) - global alignment function
""" shuffled_seqs = []
for iterations in range(1, 100):
    shuffled_seq = np.random.permutation(np.array(list(str(mystery_seq.seq))))
    shuffled_seq = Seq(''.join(shuffled_seq))
    shuffled_seqs.append(shuffled_seq) """

# perform alignment for the mystery sequence against all sequences in the database
for i in dog_database.index:
    print("Now checking alignment with:", i, dog_database.loc[i, "breed"]) # progress bar
    alignment = aligner.align(dog_database.loc[i, "sequence"], mystery_seq.seq) # Needleman-Wunsch/global alignment
    dog_database.loc[i, "align_score"] = alignment.score # writes alignment score to dataframe, saves
    
    # Local alignment only: Karlin-Altschul Algorithm (E- and P-value for statistical quality of alignment, same algorithm as BLAST)
    score = alignment.score 
    e_value = k_value * length_seq * len(dog_database.index) * math.exp(-l * score) # probability of getting a score more than the one you generate above
    p_value = 1 - math.exp(-e_value) # estimates p-value from e-value. Assumes a Gumbel extreme value distribution of sequences.
    print(e_value, p_value)
    dog_database.loc[i, "E_value"] = e_value # writes alignment score to dataframe, saves

    # calculate alignment scores for the shuffles mystery sequences, find p-value based on this - global alignment function
    # (p-value = proportion of randomly shuffled sequence scores being equal or greater than the original alignment score)
    # NOTE: Very slow, find a more efficient way to generate p-values?
"""     for j in shuffled_seqs:
        shuffled_scores = []
        shuffled_alignment = aligner.align(dog_database.loc[i, "sequence"], j)
        shuffled_scores.append(shuffled_alignment.score)
        p_value = sum(score >= dog_database.loc[i, "align_score"] for score in shuffled_scores) / 100
        dog_database.loc[i, "p_value_randomseq"] = p_value """

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