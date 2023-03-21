import similarity_checker as sim

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

# def alignment
def test_alignment():
    seq2 = Seq("AATTCCCCGG") # 100% match with 1
    seq3 = Seq("GGGGGGGGCC") # Low match with 1, 2, 4
    seq4 = Seq("GATTCCCCGG") # High match with 1 and 2, low match with 4
    pass