from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Align
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo import TreeConstruction as Tree
from Bio.Phylo.Consensus import *
import math
import os
import pandas as pd
import matplotlib.pyplot as plt

# To run: cd tests, pytest .\TEST.py
# To do: Put into a test class

# Database import/global class attributes

# def __init__
def test_init():
    pass

# def alignment
def test_alignment():
    pass

seq1 = Seq("AATTCCCCGG")
seq2 = Seq("AATTCCCCGG")

aligner = Align.PairwiseAligner()
aligner.mode = "global"

alignment = aligner.align(seq1, seq2)
print(alignment.score)

# def karlin_altschul
def test_karlinaltschul():
    pass

# def phylogeny
def test_phylogeny():
    pass