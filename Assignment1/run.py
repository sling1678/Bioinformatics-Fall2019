from file_readers import read_fasta_file
from translate import translate
from secondary_structure.hydrophobicity_character import compute_hydrophobicity_character
#from secondary_structure import hydrophobicity_character
import pickle

import matplotlib.pyplot as plt

data_file = "data/Assignment1Sequences.txt"
mRNAs = read_fasta_file.read_file(data_file)
aminoacid_sequences = []
for mRNA in mRNAs:
  label, header = mRNA[0], mRNA[1]
  aminoacid_sequence, nucleotide_errors = \
    translate.translate_simple(mRNA[2], one_letter=True)
  aminoacid_sequences.append([label, header, aminoacid_sequence, nucleotide_errors])

# TODO
for index, sample in enumerate(aminoacid_sequences):
  aa_sequence = sample[2]
  hb = compute_hydrophobicity_character(aa_sequence, aa_symbol_size=1, span_size=13)
  aminoacid_sequences[index] = [sample[0], sample[1], sample[2], sample[3], hb]

# TODO - generate plots and criteria for selecting transmembrane regions 
# plt.plot(hb, "-o")
# plt.show()

# with open('data/aa_sequences.dat', 'wb') as outfile:
#     pickle.dump(aminoacid_sequences, outfile)

