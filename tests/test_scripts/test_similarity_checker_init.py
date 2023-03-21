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

class TestInit():
    seq1 = Seq("AATTCCCCGG")
    instance1 = sim.MysterySequence(seq1)

    def test_instance_attr_1(self, testinstance = instance1):
        """Checks that the class instance has a string 'self.sequence_path' attribute"""
        assert isinstance(testinstance.sequence_path, str)

    def test_instance_attr_2(self):
        pass
