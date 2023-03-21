import similarity_checker as sim

from Bio import Align
import os
import pandas as pd

class TestAlignment:
    test_seq = sim.MysterySequence("./tests/test_seq")

    seq2 = "AATTCCCCGG" # 100% match with 1
    seq3 = "GGGGGGGGCC" # Low match with 1
    seq4 = "GATTCCCCGG" # High/close match with 1
    seq5 = "XXXXXXXXXX" # No match with 1
    test_ref = pd.DataFrame(
        {
            "breed": ["seq2", "seq3", "seq4", "seq5"],
            "sequence": [seq2, seq3, seq4, seq5],
        }
    )

    refseqalign = sim.MysterySequence.alignment(test_seq, test_ref, "./tests/test_reference/")

    def test_alignment(self, sequence = test_seq, database = test_ref):
        """Checks that the code produces pairwise alignments"""
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        for i in database.index:
            alignment = aligner.align(database.loc[i, "sequence"], sequence.seq) # instance already confirmed to have .seq attribute
        assert isinstance(alignment, Align.PairwiseAlignments)

    def test_score_1(self, sequence = test_seq, database = test_ref):
        """Checks that the alignments in the function have a score"""
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        for i in database.index:
            alignment = aligner.align(database.loc[i, "sequence"], sequence.seq) # instance already confirmed to have .seq attribute
        assert isinstance(alignment.score, float)

    def test_score_2(self, output=refseqalign):
        """seq and seq2 produces an alignment score of 10 (10/10 complete alignment)"""
        assert output.loc[output["breed"] == "seq2", "raw_alignment_score"].item() == 10

    def test_score_3(self, output=refseqalign):
        """seq and seq3 produces a score below that for seq2"""
        assert output.loc[output["breed"] == "seq3", "raw_alignment_score"].item() < 10

    def test_score_4(self, output=refseqalign):
        """seq and seq4 produces a score higher than that for seq3"""
        assert output.loc[output["breed"] == "seq4", "raw_alignment_score"].item() > output.loc[output["breed"] == "seq3", "raw_alignment_score"].item()

    def test_score_5(self, output=refseqalign):
        """seq and seq5 produces an alignment score of 0 (no letters align)"""
        assert output.loc[output["breed"] == "seq5", "raw_alignment_score"].item() == 0

    def test_output_1(self, output=refseqalign):
        """The function outputs a pandas dataframe"""
        assert isinstance(output, pd.DataFrame)

    def test_output_2(self, database=test_ref, output=refseqalign):
        """The database has the right number of sequences (same as the input reference dataframe)"""
        assert output.shape[0] == database.shape[0]

    def test_output_3(self, output=refseqalign):
        """The database has the right number of columns (+1 compared to input)"""
        assert output.shape[1] == 2 + 1

    def test_output_4(self, output=refseqalign):
        """The added column contains numerical scores/data"""
        assert output["raw_alignment_score"].dtypes == float

    def test_output_5(self, output=refseqalign):
        """The dataframe has been sorted so the highest score is the top row"""
        assert output.loc[0, "raw_alignment_score"] == max(output["raw_alignment_score"])

    def test_output_6(self, output=refseqalign):
        """The printed reference sequence is the one with the highest alignment score"""
        assert output.iloc[0, 0] == "seq2"

    def test_output_7(self):
        """The dataframe is saved as a .csv to the specified directory"""
        assert os.path.exists("./tests/test_reference/similarity_alignment_raw_scores.csv")