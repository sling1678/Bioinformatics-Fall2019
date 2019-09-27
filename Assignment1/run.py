from file_readers import read_fasta_file
from translate import translate
from hydrophobicity import trapezoid_rule_based_profile
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


# Translate each gene sequence
for mRNA in mRNAs:
  # Save name of gene
  label, header = mRNA[0], mRNA[1]
  # Save translated amino acid sequence and nucleotide errors
  aminoacid_sequence, nucleotide_errors = \
    translate.translate_simple(mRNA[2])
  aminoacid_sequences.append([label, header, aminoacid_sequence, nucleotide_errors])


# Create hydrophobicity graphs
for idx, sample in enumerate(aminoacid_sequences):
  # Calculate hydrophobicity
  aa_sequence = sample[2]

  hb = trapezoid_rule_based_profile.build_hydrophobicity_profile(aa_sequence)
  aminoacid_sequences[idx].append(hb)

  # Plot the data
  xdata = np.array(range(len(hb)))
  ydata = np.array(hb)
  fig = plt.figure(figsize=(10,6), num=idx)
  ax = fig.add_subplot(1, 1, 1)
  ax.plot(xdata, ydata, '-b')
  ax.set_xlabel("Amino acid Sequence Number")
  ax.set_ylabel(f"Sum of Product of Relative Hydrophobicity and Relative Weight")
  ax.set_title(f"Hydrophobicity Plot of {sample[0]}")

  # Hydrophobicity parameter cut off (Above: 'M', Shaded: 'P', Bellow: 'x')
  y1 = 0.5    # lower cut off
  y2 = 1.0    # upper cut off
  color = 'red'    # color of shaded putative region
  ax.axhspan(y1, y2, facecolor=color, alpha=0.5)

  # Display graphs
  plt.show()


# Analyze hydrophobicity profile and save predictions of transmembrane domains
for prediction, sample in enumerate(aminoacid_sequences):
  aa_sequence = sample[2]

  hb = trapezoid_rule_based_profile.analyze_sequence(aa_sequence)
  aminoacid_sequences[prediction].append(hb)


# Output gene sequence name, translated amino acid sequences (one letter), nucleotide errors, and predicted transmembrane domains
for aminoacid_sequence in aminoacid_sequences: 
  print(aminoacid_sequence[0])
  print("\nAmino acid Sequence (one letter translation):")
  print(aminoacid_sequence[2])
  print("\nNucleotide Errors:")
  print(aminoacid_sequence[3])
  print("\nPrediction of Transmembrane Domains ('x': undecided, 'M': definitely inside membrane, 'P': putatively inside membrane):")
  print(aminoacid_sequence[5])
  print("\n")
