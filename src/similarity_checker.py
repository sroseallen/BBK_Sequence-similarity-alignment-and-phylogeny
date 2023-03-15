from Bio import SeqIO
from Bio import Align
import math
import os
import pandas as pd
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
aligner.mode = "local" # set to local to allow E-value calculation (model for E-value is based on local alignments)

# set up Karlin-Altschul parameters for E-value and P-value generation
k = 0.1 # default value for k, normalising constant. Given in 'BLAST: An essential guide to the Basic Local Alignment Tool, 2003'
m = len(mystery_seq)
all_database_bases = "".join(dog_database["sequence"])
n = len(all_database_bases)

# frequency of all bases in the mystery sequence
mystery_freq = {"A": mystery_seq.count("A") / m,
                "C": mystery_seq.count("C") / m,
                "T": mystery_seq.count("T") / m,
                "G": mystery_seq.count("G") / m}
p1 = sum(mystery_freq.values())

# estimation of lambda value given frequency of bases in all of dog_database
dog_freq = {"A": all_database_bases.count("A") / n,
            "C": all_database_bases.count("C") / n,
            "T": all_database_bases.count("T") / n,
            "G": all_database_bases.count("G") / n}
p2 = sum(dog_freq.values())

# estimated score if predicted a perfect match - ie the % of bases in mystery_seq which would on average in the database align exactly like A:A, G:G, C:C or T:T
A_A = mystery_freq["A"] * dog_freq["A"]
C_C = mystery_freq["C"] * dog_freq["C"]
G_G = mystery_freq["G"] * dog_freq["G"]
T_T = mystery_freq["T"] * dog_freq["T"]

average_predicted = A_A + C_C + G_G + T_T

lambda_ka = math.log(1 / p1*p2) / average_predicted # equation to get Sλ, normalised score in K-A equation (converted from Sλ to λ alone by dividing by combined predicted frequency of perfectly aligned bases)

# perform alignment for the mystery sequence against all sequences in the database
for i in dog_database.index:
    print("Now checking alignment with:", i, dog_database.loc[i, "breed"]) # progress bar
    alignment = aligner.align(dog_database.loc[i, "sequence"], mystery_seq.seq) # Smith-Waterman local alignment
    dog_database.loc[i, "raw_alignment_score"] = alignment.score # writes alignment score to dataframe, saves

    # generate E and P values
    e_value = k * m * n * math.exp(-(lambda_ka * alignment.score)) # probability of getting a score higher than the one you generate above with a random sequence
    p_value = 1 - math.exp(-e_value) # estimates p-value from e-value. Assumes a Gumbel extreme value distribution of sequences.
    dog_database.loc[i, "E_value"] = e_value # writes score to dataframe
    dog_database.loc[i, "p_value"] = p_value # writes score to dataframe

# output table of aligned scores sorted by alignment score
dog_database = dog_database.sort_values("p_value", ascending=False)
dog_database.to_csv("results/final/similarity_alignment.csv")

# phylogeny tree
# global multiple sequence alignment for database of sequences
# calculator = DistanceCalculator('identity') - identity is the scoring model, but can check attr of calculator to see other models available
# above.get_distance(aln) - gives a DistanceMatrix object
# constructor = DistanceTreeConstructor(calculator, 'nj') - uses nearest neighbour joining method
# tree = constructor.build_tree(aln)
# bootstrapping to get the best tree: 
# trees = bootstrap_trees(aln, 100, constructor, majority_consensus) - generates just the consensus tree
# draw() to get a matplotlib rooted tree image output