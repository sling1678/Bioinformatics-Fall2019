from file_readers import read_fasta_file
from translate import translate
# from secondary_structure.hydrophobicity_character import compute_hydrophobicity_character
from secondary_structure import hydrophobicity_character
import pickle
import numpy as np 

import matplotlib.pyplot as plt

data_file = "data/Assignment1Sequences.txt"
mRNAs = read_fasta_file.read_file(data_file)
aminoacid_sequences = []
for mRNA in mRNAs:
  label, header = mRNA[0], mRNA[1]
  # Just to check length of genome
  aminoacid_sequence, nucleotide_errors = \
    translate.translate_simple(mRNA[2], one_letter=True)
  aminoacid_sequences.append([label, header, aminoacid_sequence, nucleotide_errors])

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


# save 
# with open('data/aa_sequences.dat', 'wb') as outfile:
#     pickle.dump(aminoacid_sequences, outfile)

