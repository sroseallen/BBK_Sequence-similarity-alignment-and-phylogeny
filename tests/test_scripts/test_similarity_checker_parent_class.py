import similarity_checker as sim

from Bio.SeqRecord import SeqRecord
import os
import pandas as pd
import pytest

class TestParentClass:
    def test_database_file(self):
        """Test 1.0: Raises an error if the reference file is missing"""
        with pytest.raises(FileNotFoundError):
            sim.MysterySequence("./tests/test_empty")

    def test_database_list_1(self):
        """Test 1.1: The correct reference database file is in the ./data file"""
        path = "./data/reference/dog_breeds.fa"
        assert os.path.exists(path)

    def test_database_list_2(self):
        """Test 1.2: The files in the data directory are definitely fasta files"""
        assert [".fa" in x for x in os.listdir("./data/reference")]

    def test_database_list_3(self):
        """Test 1.3: The database has been successfully imported a SeqRecord object"""
        database_list = sim.MysterySequence.database_list
        assert isinstance(database_list[0], SeqRecord)

    def test_description_list_1(self):
        """Test 1.4: The database converts to a list in the parsing loop"""
        description_list = sim.MysterySequence.description_list
        assert isinstance(description_list, list)

    def test_description_list_2(self):
        """Test 1.5: The database description contains no '[' characters"""
        description_list = sim.MysterySequence.description_list
        assert "[" not in description_list

    def test_breed_1(self):
        """Test 1.6: There is a breed listed for every sequence in the reference database"""
        reference_breeds = sim.MysterySequence.reference_breeds
        assert len(reference_breeds) == len(sim.MysterySequence.database_list)

    def test_breed_2(self):
        """Test 1.7: The breeds being added are character strings"""
        reference_breeds = sim.MysterySequence.reference_breeds
        assert [isinstance(x, str) for x in reference_breeds]

    def test_breed_3(self):
        """Test 1.8: The breeds being added are actually breeds (not sample IDs)"""
        reference_breeds = sim.MysterySequence.reference_breeds
        assert ["|" not in x for x in reference_breeds]

    def test_sample_1(self):
        """Test 1.9: There is a breed listed for every sequence in the reference database"""
        reference_sample_id = sim.MysterySequence.reference_sample_id
        assert len(reference_sample_id) == len(sim.MysterySequence.database_list)

    def test_sample_2(self):
        """Test 1.10: The sample IDs are actually sample IDs"""
        reference_sample_id = sim.MysterySequence.reference_sample_id
        assert ["gb" in x for x in reference_sample_id]

    def test_sequence_1(self):
        """Test 1.11: There is a sequence listed for every sequence in the reference database"""
        reference_sequence = sim.MysterySequence.reference_sequence
        assert len(reference_sequence) == len(sim.MysterySequence.database_list)

    def test_sequence_2(self):
        """Test 1.12: The sequences are parsed as strings"""
        reference_sequence = sim.MysterySequence.reference_sequence
        assert [isinstance(x, str) for x in reference_sequence]

    def test_sequence_3(self):
        """Test 1.13: The sequences are gapped nucleotide mitochondrial sequences"""
        reference_sequence = sim.MysterySequence.reference_sequence
        for x in reference_sequence:
            for letter in x:
                assert letter in 'ACTG-NR' # -, N, and R are artefacts in the reference database - expected to be there.

    def test_reference_database_1(self):
        """Test 1.14: The reference db is a pandas dataframe"""
        reference_database = sim.MysterySequence.reference_database
        assert isinstance(reference_database, pd.DataFrame)

    def test_reference_database_2(self):
        """Test 1.15: The dataframe has 4 columns, all the length of the number of sequences in the original file"""
        reference_database = sim.MysterySequence.reference_database
        assert reference_database.shape == (len(sim.MysterySequence.database_list), 4)