from file_readers import read_fasta_file
from translate import translate
# from secondary_structure.hydrophobicity_character import compute_hydrophobicity_character
from secondary_structure import hydrophobicity_character
import pickle
import numpy as np
import matplotlib.pyplot as plt

"""
BIOINFORMATICS - Assignment 1

Group: David Blanck, Sammy Ling, Andrew Richter, Kaveri Sharma

Description: Program for taking gene sequences from a FASTA file and translating into amino acid sequence.
  Outputs translated sequences.
  Computes hydrophobic average over span.
  Will determine how likely each sequence is to encode a membrane protein.

"""

data_file = "data/Assignment1Sequences.txt"    # sample gene sequence file
mRNAs = read_fasta_file.read_file(data_file)    # read in file and store genes to be translated
aminoacid_sequences = []    # translated into one letter amino acid sequences
pretty_aminoacid_sequences = []    # translated in three letter amino acid sequences


# translate each gene
for mRNA in mRNAs:
  label, header = mRNA[0], mRNA[1]
  # Just to check length of genome
  aminoacid_sequence, nucleotide_errors = \
    translate.translate_simple(mRNA[2], one_letter=True)
  aminoacid_sequences.append([label, header, aminoacid_sequence, nucleotide_errors])

  # generate three letter sequence for output
  pretty_aminoacid_sequence, nucleotide_errors = \
    translate.translate_simple(mRNA[2], one_letter=False)
  pretty_aminoacid_sequences.append([label, header, pretty_aminoacid_sequence, nucleotide_errors])


# create hydrophobicity graphs
span_size=13
for idx, sample in enumerate(aminoacid_sequences):
  aa_sequence = sample[2]

  hb = hydrophobicity_character.compute_hydrophobicity_character(aa_sequence, aa_symbol_size=1, span_size=span_size)
  aminoacid_sequences[idx].append(hb)

  # plot the data
  xdata = np.array(range(len(hb))) + span_size//2
  ydata = np.array(hb)
  fig = plt.figure(figsize=(10,6), num=idx)
  ax = fig.add_subplot(1, 1, 1)
  ax.plot(xdata, ydata, '-b')
  ax.set_xlabel("Aminoacid Sequence Number")
  ax.set_ylabel(f"Hydrophobicity ove span of {span_size} aminoacids")
  ax.set_title(f"Hydrophobicity Plot of {sample[0]}")
  plt.show()


# Output gene sequence, gene length, aa lengthtranslated amino acid sequences (one & three letter), and nucleotide errors
for (mRNA, aminoacid_sequence, pretty_aminoacid_sequence) in zip(mRNAs, aminoacid_sequences, pretty_aminoacid_sequences): 
  print(aminoacid_sequence[0])
  print("\nGene Sequence Length:")
  print(len(mRNA[2]))
  print("\nAmino acid Sequence Length:")
  print(len(aminoacid_sequence[2]))
  print("\nAmino acid Sequence (three letter translation):")
  print(pretty_aminoacid_sequence[2])
  print("\nAmino acid Sequence (one letter translation):")
  print(aminoacid_sequence[2])
  print("\nNucleotide Errors:")
  print(aminoacid_sequence[3])
  print("\n")


# save
# with open('data/aa_sequences.dat', 'wb') as outfile:
#     pickle.dump(aminoacid_sequences, outfile)
