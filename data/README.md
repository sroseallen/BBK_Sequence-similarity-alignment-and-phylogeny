# data: database, and test sequence

The `dog_breeds.fa` file contains the sequence in our database, to compare against.
The `mystery.fa` sequence contains a sequence from an uknown dog type, and we want to identify what is the closest breed to it in our database.


# Phylogeny

If you want to make a phylogeny, the simplest way to do it is to feed the file to an online webservice, such as https://www.ebi.ac.uk/Tools/phylogeny/simple_phylogeny/
export a file, then read it in using the Phylo submodule.

You can also compute a phylogenetic tree via the following two steps:
- construct a distance matrix
- construct a tree using a construction algorithm
See the following page for examples: https://biopython.org/docs/1.75/api/Bio.Phylo.TreeConstruction.html

