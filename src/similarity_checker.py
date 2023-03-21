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


class MysterySequence:
    """
    A class to handle the unknown sequence to be aligned and described against the reference database.
    Contains class attributes/methods for the reference dog database set.

    Class Attributes:
        database_list(List): Reads in all data from the reference database in the ./data file of this module.
        reference_breeds(List): List of breeds in the reference database
        reference_sample_id(List) List of samples in the reference database
        reference_sequence(List): List of all sequences in the reference database
        reference_database(pd.DataFrame): A pandas dataframe of all sequences in the reference database (includes breed name, sample id, and sequence)

    Object Attributes:
        seq (Bio.Seq): Sequence object for the unknown sequence
        sequence_path(str): Filepath for the unknown sequence

    Object Methods:
        __init__(self, seq, sequence_directory="./data/seq"): Create a MysterySequence object
        alignment(self, database=reference_database, save_dir="results/final/"): Align unknown sequence to all sequences in the reference database
        karlin_altschul(self, aligned_database, k:float=0.1, save_dir="results/final/"): Add E- and P- values to all aligned raw scores
        phylogeny(self:, database=reference_database, bootstrap_iteration=1, save_dir="results/final/"): Produce a consensus phylogeny tree for the ref database and unknown sequence
    """
    try:
        database_list = list(
            SeqIO.parse("./data/reference/dog_breeds.fa", "fasta")
        )
    except:
        raise FileNotFoundError ("No reference file found - please ensure there is a .fasta format dog_breeds database file in the ./data/reference directory.")
    
    reference_breeds = []  # initialise lists for input into pandas dataframe
    reference_sample_id = []
    reference_sequence = []

    for (
        dog
    ) in (
        database_list
    ):
        description_list = dog.description.split("[")
        breed = description_list[7][6:-2]
        reference_breeds.append(
            breed
        )  # adds data for each dog to the dataframe row by row
        reference_sample_id.append(dog.id)
        reference_sequence.append("".join(dog.seq))

    reference_database = pd.DataFrame(
        {
            "breed": reference_breeds,
            "sample_id": reference_sample_id,
            "sequence": reference_sequence,
            "align_score": 0,
        }
    )

    def __init__(self, sequence_path: str = "./data/seq") -> "MysterySequence":
        """
        Creates a MysterySequence instance from a .fasta file (path to the .fasta file defined in sequence_directory)
        """
        self.sequence_path = sequence_path
        if len(os.listdir(sequence_path)) > 1:
            raise OSError ("Too many files! Please ensure there is only one fasta file in the sequence directory")
        if len(os.listdir(sequence_path)) == 0:
            raise FileNotFoundError ("No file found! Please ensure your .fasta sequence file is saved in the specified directory.")
        for filename in os.listdir(
            sequence_path
        ):  # generic to allow any named fasta file
            path = os.path.join(sequence_path, filename)
            if filename.endswith(".fa") or filename.endswith(".fasta"):
                read = SeqIO.read(path, "fasta")
                self.seq = read.seq
            else:
                raise FileNotFoundError ("No fasta file found, please ensure your .fasta sequence file is saved in the specified directory.")

    def alignment(
        self: "MysterySequence",
        database: pd.DataFrame = reference_database,
        save_dir: str = "./results/final/",
    ) -> pd.DataFrame:
        """
        Performs a pairwise alignment of the unknown sequence to each sequence in the class-defined reference database.

        Arguments:
            self: An object of the class MysterySequence
            database: A pandas dataframe defined for all class objects, defaults to the dog sequence reference database for this class
            save_dir: Path to the directory/folder the output should be saved in (defaults to the "results/final/" folder in the package directory)
        """
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"

        # perform alignment for the mystery sequence against all sequences in the database
        for i in database.index:
            print(
                "Now checking alignment against:", i, database.loc[i, "breed"]
            )  # progress bar
            alignment = aligner.align(
                database.loc[i, "sequence"], self.seq
            )  # Global alignment
            database.loc[
                i, "raw_alignment_score"
            ] = alignment.score  # writes alignment score to dataframe

        # output table of aligned scores sorted by alignment score
        aligned_database = database.sort_values("raw_alignment_score", ascending=False)
        print("The most similar sequence is", aligned_database.iloc[0, 0])
        save_path = save_dir + "similarity_alignment_raw_scores.csv"
        aligned_database.to_csv(save_path)

        return aligned_database

    def karlin_altschul(
        self: "MysterySequence",
        aligned_database: pd.DataFrame,
        k: float = 0.1,
        save_dir: str = "./results/final/",
    ) -> pd.DataFrame:
        """
        For a saved aligned file, performs Karlin-Altschul statistical test for each alignment, and records E-value and P-value in the database.

        Arguments:
            self: An object of the class MysterySequence
            aligned_database: A pandas dataframe, which MUST include at minimum the reference sequences and the raw alignment scores for each sequence.
            k: Normalising constant for the Karlin-Altschul test. Defaults to 0.1, value provided in 'BLAST: An essential guide to the Basic Local Alignment Tool, 2003'
            save_dir: Path to the directory/folder the output should be saved in (defaults to the "results/final/" folder in the package directory)
        """
        m = len(self.seq)
        all_database_bases = "".join(aligned_database["sequence"])
        n = len(all_database_bases)

        # frequency of all bases in the mystery sequence
        mystery_freq = {
            "A": self.seq.count("A") / m,
            "C": self.seq.count("C") / m,
            "T": self.seq.count("T") / m,
            "G": self.seq.count("G") / m,
        }
        p1 = sum(mystery_freq.values())

        # frequency of all bases in the reference database
        dog_freq = {
            "A": all_database_bases.count("A") / n,
            "C": all_database_bases.count("C") / n,
            "T": all_database_bases.count("T") / n,
            "G": all_database_bases.count("G") / n,
        }
        p2 = sum(dog_freq.values())

        # estimated score if predicted a perfect match - ie the % of bases in mystery_seq which would on average in the database align exactly like A:A, G:G, C:C or T:T
        A_A = mystery_freq["A"] * dog_freq["A"]
        C_C = mystery_freq["C"] * dog_freq["C"]
        G_G = mystery_freq["G"] * dog_freq["G"]
        T_T = mystery_freq["T"] * dog_freq["T"]

        average_predicted = A_A + C_C + G_G + T_T

        lambda_ka = (
            math.log(1 / p1 * p2) / average_predicted
        )  # equation to get Sλ, normalised score in K-A equation (converted from Sλ to λ alone by dividing by combined predicted frequency of perfectly aligned bases)

        # generate E (Expect) value and P value
        for i in aligned_database.index:
            e_value = (
                k
                * m
                * n
                * math.exp(
                    -(lambda_ka * aligned_database.loc[i, "raw_alignment_score"])
                )
            )  # probability of getting a score higher than the one you generate above with a random sequence
            p_value = 1 - math.exp(
                -e_value
            )  # estimates p-value from e-value. Assumes a Gumbel extreme value distribution of sequences.
            aligned_database.loc[i, "E_value"] = e_value  # writes score to dataframe
            aligned_database.loc[i, "p_value"] = p_value  # writes score to dataframe

        aligned_database = aligned_database.sort_values("E_value", ascending=True)
        print("The most statistically likely alignment is against", aligned_database.iloc[0, 0])
        save_path = save_dir + "similarity_alignment_statistical_test.csv"
        aligned_database.to_csv(save_path)

        return aligned_database

    def phylogeny(
        self: "MysterySequence",
        database: pd.DataFrame = reference_database,
        bootstrap_iteration: int = 1,
        save_dir: str = "./results/final/",
    ) -> Phylo:
        """
        Produces a phylogenetic tree including the entire reference database and the unknown sequence.

        Arguments:
            self: An object of the class MysterySequence
            database: A pandas dataframe defined for all class objects (defaults to the dog sequence reference database for this class)
            bootstrap_iteration: Number of repeat trees created in generating the consensus tree (defaults to 1)
            save_dir: Path to the directory/folder the outputs should be saved in (defaults to the "results/final" folder in the package directory)
        """
        # global multiple sequence alignment for database of sequences
        seqs = []
        seqs.append(SeqRecord(self.seq, id="UNKNOWN SEQUENCE"))
        for i in database.index:
            info = (database.loc[i, "sample_id"]) + (database.loc[i, "breed"])
            a = SeqRecord(Seq(database.loc[i, "sequence"]), id=info)
            seqs.append(a)
        multi_align = Align.MultipleSeqAlignment(seqs)

        save_path = save_dir + "database_mulitpleseq_alignment.clustal"
        AlignIO.write(multi_align, save_path, "clustal")

        # Distance Calculator - identity is the scoring model (default, works for both DNA and protein sequences)
        calculator = Tree.DistanceCalculator("identity")

        # bootstrapping to get the best tree:
        # "trees" can produce replicate trees if bootstrapping_iteration is provided, and returns just the main consensus tree (uses nearest neighbour joining method)
        boot_constructor = Tree.DistanceTreeConstructor(calculator, "nj")
        consensus_tree = bootstrap_consensus(
            multi_align, bootstrap_iteration, boot_constructor, majority_consensus
        )

        # draw() to get a matplotlib rooted tree image output (networkx not used as BioPhylo docs suggest branch lengths not indicative of evolutionary distance)
        consensus_tree.rooted = True
        consensus_tree.root.color = "green"

        # write the consensus tree in multiple formats
        save_path = save_dir + "consensus_tree.xml"
        Phylo.write(consensus_tree, save_path, "phyloxml")

        fig = plt.figure(figsize=(10, 20), dpi=600)
        axes = fig.add_subplot(1, 1, 1)
        Phylo.draw(consensus_tree, axes=axes, do_show=False)
        save_path = save_dir + "consensus_tree.png"
        plt.savefig(save_path)

        return


if __name__ == "__main__":
    unknown = MysterySequence()
    output1 = unknown.alignment()
    output2 = unknown.karlin_altschul(aligned_database=output1)
    output3 = unknown.phylogeny()
