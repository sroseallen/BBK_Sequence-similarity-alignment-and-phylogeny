from Bio import SeqIO
from Bio import Align
import os
import scipy as sc
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
# tests: length of sequence = the same? 
alignment_scoring = []
for i in dog_database.index:
    print("Now checking alignment with:", i, dog_database.loc[i, "breed"]) # progress bar
    aligner = Align.PairwiseAligner()
    alignment = aligner.align(dog_database.loc[i, "sequence"], mystery_seq.seq) # Needleman-Wunsch/global alignment
    dog_database.loc[i, "align_score"] = alignment.score # writes alignment score to dataframe, saves
dog_database = dog_database.sort_values("align_score")

dog_database.to_csv("results/final/similarity_alignment.csv")

# probability scores
# scipy.stats import pearsonr (standard correlation coefficient)
# def pearsonr_pval(x,y):return pearsonr(x,y)[1]
# corr = df.corr(method=pearsonr_pval) - generates a correlation matrix for all sequences vs all sequences in the database
# perhaps more usefully: define col in dataframe for p-value
# then, for each row, call pearsonr_pval on the alignment score against the maximum possible alignment score (ie: alignment of the mystery sequence with itself)
# hopefully, this would then return a pearson correlation p-value for that sequence against the mystery sequence?

# phylogeny tree
# global multiple sequence alignment for database of sequences
# calculator = DistanceCalculator('identity') - identity is the scoring model, but can check attr of calculator to see other models available
# above.get_distance(aln) - gives a DistanceMatrix object
# constructor = DistanceTreeConstructor(calculator, 'nj') - uses nearest neighbour joining method
# tree = constructor.build_tree(aln)
# bootstrapping to get the best tree: 
# trees = bootstrap_trees(aln, 100, constructor, majority_consensus) - generates just the consensus tree
# draw() to get a matplotlib rooted tree image output