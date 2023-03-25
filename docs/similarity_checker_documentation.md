# About Similarity Checker
This module is based around the 'MysterySequence' class, which can be initialised for any fasta file containing a single sequence. Information about the class can be found below, or by using the help() function in your python environment.

There are two main ways to use this module:
1. Running the whole script directly from the command line. This uses all the default values, which have been set to provide all outputs in-line with the Biocomputing coursework (eg: will save all results in the 'results' directory, will use the unknown sequence and reference database saved under the 'data' directory). 

2. Import of the module and running each function separately. This allows specification of non-default values and for sequences other than the unknown sequence provided for the Biomputing course to be run with the module. As all functions are actually methods of the MysterySequence class, the unknown sequence must first be instantiated using 'instance = MysteryChecker()'.

# Class: similarity_checker.MysterySequence()
A class to handle the unknown sequence to be aligned and described against the reference database. Every instance of this class can use the reference dog sequence database in its methods, allowing for easy alignment and calling of values in the reference database. 

    Class Attributes:
        database_list(List): Reads in all data from the reference database in the .\data file of this module.
        reference_breeds(List): List of breeds in the reference database
        reference_sample_id(List) List of samples in the reference database
        reference_sequence(List): List of all sequences in the reference database
        reference_database(pd.DataFrame): A pandas dataframe of all sequences in the reference database (includes breed name, sample id, and sequence)

    Object Attributes:
        seq (Bio.Seq): Sequence object for the unknown sequence
        sequence_path(str): Filepath for the unknown sequence
    
# Functions / Methods
There are four methods for the MysterySequence class, which are detailed below (and also by calling help() for any function)

## __init__ (self, seq, sequence_directory=".\data\seq")
Creates a MysterySequence instance from a .fasta file (path to the .fasta file defined in sequence_directory). Must be a .fasta file, otherwise the class instance will not be created.

Arguments:
  *self: Instance name
  *seq: Instance attribute containing the sequence information from the .fasta file
  *sequence_directory: Where the unknown sequence is saved - this defaults to the data folder in the project directories, but can be specified to any path accessible.

## .alignment (self, database=reference_database, save_dir="results/final/")
Performs a pairwise alignment of the unknown sequence to each sequence in the class-defined reference database.

Arguments:
  *self: An object of the class MysterySequence
  *database: A pandas dataframe defined for all class objects, defaults to the dog sequence reference database for this class
  *save_dir: Path to the directory/folder the output should be saved in (defaults to the "results/final/" folder in the package directory)

## .karlin_altchul (self, aligned_database, k:float=0.1, save_dir="results/final/")
For a saved aligned file, performs Karlin-Altschul statistical test for each alignment, and records E-value and P-value in the database.

Arguments:
  *self: An object of the class MysterySequence
  *aligned_database: A pandas dataframe, which MUST include at minimum the reference sequences and the raw alignment scores for each sequence. 
  *k: Normalising constant for the Karlin-Altschul test. Defaults to 0.1, value provided in 'BLAST: An essential guide to the Basic Local Alignment Tool, 2003'
  *save_dir: Path to the directory/folder the output should be saved in (defaults to the "results/final/" folder in the package directory)

## .phylogeny (self:, database=reference_database, bootstrap_iteration=1, save_dir="results/final/")
Produces a phylogenetic tree including the entire reference database and the unknown sequence.

Arguments:
  *self: An object of the class MysterySequence
  *database: A pandas dataframe defined for all class objects (defaults to the dog sequence reference database for this class)
  *bootstrap_iteration: Number of repeat trees created in generating the consensus tree (defaults to 1 - no bootstrapping enabled by default).
  *save_dir: Path to the directory/folder the outputs should be saved in (defaults to the "results/final" folder in the package directory)