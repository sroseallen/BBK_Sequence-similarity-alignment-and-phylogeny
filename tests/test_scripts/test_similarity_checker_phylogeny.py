import similarity_checker as sim

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Align
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo import TreeConstruction as Tree
from Bio.Phylo.Consensus import *
import os
import pandas as pd
import matplotlib.pyplot as plt

class TestPhylogeny:
    """
    Test class for the similarity_checker.MysterySequence.phylogenyfunction.

    Class Attributes:
        test_seq, seq2-5: Sequences to make up a mystery test sequence and 4 further reference sequences, each is 10 characters long.
        test_ref: A test reference database made of 4 sequences of 10 characters each (from seq2-seq5).
        phyl: Output from running the similarity_checker.MysterySequence.phylogeny function with the test attributes above

    Object Methods (tests run):
        test_multiseq_1: Mystery Sequence can convert to a SeqRecord instance
        test_multiseq_2: Info for each sequence contains both the sample_id and breed (seq2 checked)
        test_multiseq_3: Produces multiple sequence alignment object
        test_output_1: The clustal alignment is saved to the specified directory
        test_output_2: The phylogeny .xml alignment is saved to the specified directory
        test_output_3: The phylogeny alignment is saved as a static .png to the specified directory
        test_seq_distance_1: The nearest neighbour to the mystery sequence is seq2 (100% alignment)
        tset_seq_distance_2: The furthest neighbour to the mystery sequence is seq5 (0% alignment)
    """

    test_seq = sim.MysterySequence("./tests/test_seq")
    seq2 = "AATTCCCCGG" # 100% match with test_seq
    seq3 = "GGGGGGGGCC" # Low match with test_seq
    seq4 = "GATTCCCCGG" # High/close match with test_seq
    seq5 = "XXXXXXXXXX" # No match with test_seq
    test_ref = pd.DataFrame(
        {
            "breed": ["seq2", "seq3", "seq4", "seq5"],
            "sample_id": ["2", "3", "4", "5"],
            "sequence": [seq2, seq3, seq4, seq5]
        }
    )
    phyl = sim.MysterySequence.phylogeny(test_seq, test_ref, save_dir="./tests/test_reference/")

    def test_multiseq_1(self, seq=test_seq.seq):
        """Mystery Sequence can convert to a SeqRecord instance"""
        assert isinstance((SeqRecord(seq, id="unknown")), SeqRecord)
    
    def test_multiseq_2(self, database = test_ref):
        """Info for each sequence contains both the sample_id and breed (seq2 checked)"""
        info = (database.loc[0, "sample_id"]) + (database.loc[0, "breed"])
        assert info == "2seq2"

    def test_multiseq_3(self, seq=test_seq.seq, database=test_ref):
        """Produces multiple sequence alignment object"""
        seqs = [SeqRecord(seq),
                SeqRecord(Seq(database.loc[0, "sequence"])),
                SeqRecord(Seq(database.loc[1, "sequence"])),
                SeqRecord(Seq(database.loc[2, "sequence"])),
                SeqRecord(Seq(database.loc[3, "sequence"]))
                ]
        align = Align.MultipleSeqAlignment(seqs)
        assert isinstance(align, Align.MultipleSeqAlignment)

    def test_output_1(self):
        """The clustal alignment is saved to the specified directory"""
        assert os.path.exists("./tests/test_reference/database_mulitpleseq_alignment.clustal")
    
    def test_output_2(self):
        """The phylogeny .xml alignment is saved to the specified directory"""
        assert os.path.exists("./tests/test_reference/consensus_tree.xml")
    
    def test_output_3(self):
        """The phylogeny alignment is saved as a static .png to the specified directory"""
        assert os.path.exists("./tests/test_reference/consensus_tree.png")

    def test_seq_distance_1(self, output=phyl):
        """The nearest neighbour to the mystery sequence is seq2 (100% alignment)"""
        pass

    def test_seq_distance_2(self, output=phyl):
        """The furthest neighbour to the mystery sequence is seq5 (0% alignment)"""
        pass