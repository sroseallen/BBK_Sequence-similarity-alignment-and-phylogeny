import similarity_checker as sim

from Bio.Seq import Seq
import os
import pytest

class TestInit:
    instance1 = sim.MysterySequence("./tests/test_seq")

    def test_seqfile_1(self):
        """The default directory only has one sequence in it"""
        assert len(os.listdir("./data/seq")) == 1

    def test_seqfile_2(self):
        """Raises an error if there are multiple files in the sequence directory"""
        with pytest.raises(OSError):
            sim.MysterySequence("./tests/test_scripts")
    
    def test_seqfile_3(self):
        """Raises an error if there are no files in the sequence directory"""
        with pytest.raises(FileNotFoundError):
            sim.MysterySequence("./tests/test_empty")
    
    def test_seqfile_4(self):
        """Raises an error if there are no fasta files specifically in the sequence directory"""
        with pytest.raises(FileNotFoundError):
            sim.MysterySequence("./tests/test_nonfasta")

    def test_instance_attr_1(self, testinstance = instance1):
        """Checks that the class instance has a string 'self.sequence_path' attribute"""
        assert isinstance(testinstance.sequence_path, str)
    
    def test_instance_attr_2(self, testinstance = instance1):
        """Checks that the instance has a 'self.seq' attribute"""
        assert isinstance(testinstance.seq, Seq)
    
    def test_instance_class(self, testinstance = instance1):
        """Checks that the instance is successfully a MysterySequence instance"""
        assert isinstance(testinstance, sim.MysterySequence)

    
