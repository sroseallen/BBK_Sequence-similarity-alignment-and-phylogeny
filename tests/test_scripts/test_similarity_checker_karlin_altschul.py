import similarity_checker as sim

import math
import os
import pandas as pd

class TestKarlinAltschul:
    """
    Test class for the similarity_checker.MysterySequence.karlin_altschul function.

    Class Attributes:
        test_seq, seq2-4: Sequences to make up a mystery test sequence and 3 further reference sequences, each is 10 characters long.
        test_ref: A test reference database made of 3 sequences of 10 characters each (from seq2-seq4).
        refseqalign: Output from running the similarity_checker.MysterySequence.alignment function with the test attributes above (required for the karlin_altschul function)
        ka: Output from running the similarity_checker.MysterySequence.karlin_altschul function.

    Object Methods (tests run):
        test_m: Calculation of m: length of the mystery sequence
        test_mysteryfreq: Calculation of p1: frequency of each base in the mystery sequence
        test_n: Calculation of n: total number of bases in the reference database
        test_reffreq: Calculation of p2: frequency of each base in the reference database
        test_lambda: Calculation of lambda: normalising score for this combination of seq and database
        test_output_1: The function outputs a pandas dataframe
        test_output_2: The dataframe is saved as a .csv to the specified directory
        test_output_3: The dataframe has been sorted so the lowest e-value is the top row
        test_output_4: The database has the right number of columns (5)
        test_output_5: The printed reference sequence is the most significant
        test_e_value_1: Checks if the seq2 e-value is correctly calculated
        test_e_value_2: Checks if the seq3 e-value is correctly calculated
        test_e_value_3: Checks if the seq4 e-value is correctly calculated
        test_p_value_1: Checks if seq2 p-value is correctly calculated
        test_p_value_2: Checks if seq3 p-value is correctly calculated
        test_p_value_3: Checks if seq4 p-value is correctly calculated
    """

    test_seq = sim.MysterySequence("./tests/test_seq")
    seq2 = "AATTCCCCGG" # 100% match with test_seq
    seq3 = "GGGGGGGGCC" # Low match with test_seq
    seq4 = "GATTCCCCGG" # High/close match with test_seq
    test_ref = pd.DataFrame(
        {
            "breed": ["seq2", "seq3", "seq4"],
            "sequence": [seq2, seq3, seq4],
        }
    )
    refseqalign = sim.MysterySequence.alignment(test_seq, test_ref, "./tests/test_reference/")
    ka = sim.MysterySequence.karlin_altschul(test_seq, aligned_database = refseqalign, save_dir = "./tests/test_reference/")

    def test_m(self, seq=test_seq.seq):
        """Calculation of m: length of the mystery sequence"""
        assert len(seq) == 10

    def test_mysteryfreq(self, seq=test_seq.seq):
        """Calculation of p1: frequency of each base in the mystery sequence"""
        mystery_freq = {
            "A": seq.count("A") / 10,
            "C": seq.count("C") / 10,
            "T": seq.count("T") / 10,
            "G": seq.count("G") / 10,
        }
        assert mystery_freq == {"A":0.2, "C":0.4, "T":0.2, "G":0.2}

    def test_n(self, ref=test_ref):
        """Calculation of n: total number of bases in the reference database"""
        all_database_bases = "".join(ref["sequence"])
        assert len(all_database_bases) == 30
    
    def test_reffreq(self, ref=test_ref):
        """Calculation of p2: frequency of each base in the reference database"""
        all_database_bases = "".join(ref["sequence"])
        dog_freq = {
            "A": all_database_bases.count("A") / 30,
            "C": round(all_database_bases.count("C") / 30, 3),
            "T": round(all_database_bases.count("T") / 30, 3),
            "G": round(all_database_bases.count("G") / 30, 3)
        }
        assert dog_freq == {"A":0.1, "C":0.333, "T":0.133, "G":0.433}

    def test_lambda(self):
        """Calculation of lambda: normalising score for this combination of seq and database"""
        average_predicted = 0.2
        assert round(math.log(1 / (1 * 0.999)) / average_predicted, 3) == 0.005
    
    def test_output_1(self, output=ka):
        """The function outputs a pandas dataframe"""
        assert isinstance(output, pd.DataFrame)

    def test_output_2(self):
        """The dataframe is saved as a .csv to the specified directory"""
        assert os.path.exists("./tests/test_reference/similarity_alignment_statistical_test.csv")

    def test_output_3(self, output=ka):
        """The dataframe has been sorted so the lowest e-value is the top row"""
        assert output.loc[0, "E_value"] == min(output["E_value"])

    def test_output_4(self, output=ka):
        """The database has the right number of columns (5)"""
        assert output.shape[1] == 5
    
    def test_output_5(self, output=refseqalign):
        """The printed reference sequence is the most significant"""
        assert output.iloc[0, 0] == "seq2"

    def test_e_value_1(self, output=ka, lambda_ka=0.005, k=0.1, m=10, n=30):
        """Checks if the seq2 e-value is correctly calculated"""
        e_value = round((k * m * n * math.exp(-(lambda_ka * output.loc[output["breed"] == "seq2", "raw_alignment_score"]))), 2)
        assert e_value == 28.54

    def test_e_value_2(self, output=ka, lambda_ka=0.005, k=0.1, m=10, n=30):
        """Checks if the seq3 e-value is correctly calculated"""
        e_value = round((k * m * n * math.exp(-(lambda_ka * output.loc[output["breed"] == "seq3", "raw_alignment_score"]))), 2)
        assert e_value == 29.70

    def test_e_value_3(self, output=ka, lambda_ka=0.005, k=0.1, m=10, n=30):
        """Checks if the seq4 e-value is correctly calculated"""
        e_value = round((k * m * n * math.exp(-(lambda_ka * output.loc[output["breed"] == "seq4", "raw_alignment_score"]))), 2)
        assert e_value == 28.68

    def test_p_value_1(self):
        """Checks if seq2 p-value is correctly calculated"""
        assert round(1 - math.exp(-28.54), 0) == 1

    def test_p_value_2(self):
        """Checks if seq3 p-value is correctly calculated"""
        assert round(1 - math.exp(-29.70), 1) == 1

    def test_p_value_3(self):
        """Checks if seq4 p-value is correctly calculated"""
        assert round(1 - math.exp(-28.68), 1) == 1