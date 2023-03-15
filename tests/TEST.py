from Bio.Seq import Seq
from Bio import Align

seq1 = Seq("AATTCCCCGG")
seq2 = Seq("AATTCCCCGG")

aligner = Align.PairwiseAligner()
aligner.mode = "local" # set to global for the global alignment function. Local alignment is slower than global.

alignment = aligner.align(seq1, seq2)
print(alignment.score)