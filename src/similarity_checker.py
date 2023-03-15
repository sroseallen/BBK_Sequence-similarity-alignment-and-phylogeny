from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align
import Bio.Align.AlignInfo as AlgnInfo
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
# set up the aligner
# tests: length of sequence = the same? 
alignment_scoring = []
aligner = Align.PairwiseAligner()
aligner.mode = "local" # set to local to allow E-value calculation (model for E-value is based on local alignments)

# set up Karlin-Altschul parameters for E-value and P-value generation
k = 0.1 # default value for k, normalising constant. Given in 'BLAST: An essential guide to the Basic Local Alignment Tool, 2003'
m = len(mystery_seq)
n = len("".join(dog_database["sequence"]))

# frequency of all bases in the mystery sequence
mystery_freq = {"A": mystery_seq.count("A") / m,
                "C": mystery_seq.count("C") / m,
                "T": mystery_seq.count("T") / m,
                "G": mystery_seq.count("G") / m}
p1 = sum(mystery_freq.values())

# perform alignment for the mystery sequence against all sequences in the database
for i in dog_database.index:
    print("Now checking alignment with:", i, dog_database.loc[i, "breed"]) # progress bar
    alignment = aligner.align(dog_database.loc[i, "sequence"], mystery_seq.seq) # Smith-Waterman local alignment
    dog_database.loc[i, "align_score"] = alignment.score # writes alignment score to dataframe, saves
    
    # For local alignments only: Karlin-Altschul Algorithm (E- and P-value for statistical quality of alignment, same algorithm as BLASTn): Scoring Matrix Generation
    # frequency of all bases in dog_database sequence
    dog_freq = {"A" : 0, "C" : 0, "T" : 0, "G" : 0}
    dog_freq = {"A": dog_database.loc[i, "sequence"].count("A") / len(dog_database.loc[i, "sequence"]),
                "C": dog_database.loc[i, "sequence"].count("C") / len(dog_database.loc[i, "sequence"]),
                "T": dog_database.loc[i, "sequence"].count("T") / len(dog_database.loc[i, "sequence"]),
                "G": dog_database.loc[i, "sequence"].count("G") / len(dog_database.loc[i, "sequence"])}
    
    p2 = sum(dog_freq.values()) # sum of all expected frequencies, should be ~1

    # count of the observed aligned pairs
    observed_freq = {"A": {"A" : 0, "C" : 0, "T" : 0, "G" : 0}, 
                     "C": {"A" : 0, "C" : 0, "T" : 0, "G" : 0}, 
                     "T": {"A" : 0, "C" : 0, "T" : 0, "G" : 0}, 
                     "G": {"A" : 0, "C" : 0, "T" : 0, "G" : 0}}

    for i, j in zip(alignment.sequences[0], alignment.sequences[1]):
        if i not in dog_freq or j not in dog_freq:
            continue
        else:
            observed_freq[i][j] += 1

    # converting counts into frequencies
    q = 0 # sum of all observed frequencies, should be ~1
    for i in observed_freq:
        for j in observed_freq.get(i):
            observed_freq[i][j] = observed_freq[i][j] / len(mystery_seq)
            q += observed_freq[i][j]

    normalised_score = 1 / (math.log(q / (p1*p2))) # equation to get Sλ, normalised score in K-A equation
    print(normalised_score / alignment.score)

    # generate E and P values
    e_value = k * m * n * math.exp(-normalised_score) # probability of getting a score more than the one you generate above
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