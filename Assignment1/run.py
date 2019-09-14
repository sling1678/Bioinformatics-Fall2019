from file_readers import read_fasta_file
from translate import translate
import pickle

data_file = "data/Assignment1Sequences.txt"
mRNAs = read_fasta_file.read_file(data_file)
aminoacid_sequences = []
for mRNA in mRNAs:
  label, header = mRNA[0], mRNA[1]
  aminoacid_sequence, nucleotide_errors = \
    translate.translate_simple(mRNA[2], one_letter=True)
  aminoacid_sequences.append([label, header, aminoacid_sequence, nucleotide_errors])

with open('data/aa_sequences.dat', 'wb') as outfile:
    pickle.dump(aminoacid_sequences, outfile)

