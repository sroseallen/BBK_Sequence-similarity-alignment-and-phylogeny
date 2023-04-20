from src.similarity_checker import *

unknown = MysterySequence()
output1 = unknown.alignment()
output2 = unknown.karlin_altschul(output1)
output3 = unknown.phylogeny()