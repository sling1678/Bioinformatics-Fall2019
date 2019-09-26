from file_readers import read_fasta_file
from translate import translate
import data
from secondary_structure import trapezoid_rule_based_profile
import pickle
import numpy as np

data_file = "Assignment1Sequences.txt"    # sample gene sequence file
mRNAs = read_fasta_file.read_file(data_file)    # read in file and store genes to be translated
aminoacid_sequences = []    # translated into one letter amino acid sequences



# Translate each gene sequence
for mRNA in mRNAs:
  label, header = mRNA[0], mRNA[1]
  # Just to check length of genome
  aminoacid_sequence, nucleotide_errors = \
    translate.translate_simple(mRNA[2])
  aminoacid_sequences.append([aminoacid_sequence])

hydrophobicity_profiles = []
for aas in aminoacid_sequences:
  hydrophobicity_profiles.append(trapezoid_rule_based_profile.analyze_sequence(aas[0]))

print(hydrophobicity_profiles[0])
print(hydrophobicity_profiles[1])
print(hydrophobicity_profiles[2])
print(len(hydrophobicity_profiles[0]))
print(len(aas[0]))
