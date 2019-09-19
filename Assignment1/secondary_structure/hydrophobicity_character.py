import numpy as np 
from secondary_structure import aa_hydrophobicity, aminoacid_symbols
#import aa_hydrophobicity, aminoacid_symbols
def convert_to_three_letter_aa(one_letter_aa):
  aa3  = aminoacid_symbols.aminoacid_symbols()
  result = "***"
  for key, val in aa3.items():
    if one_letter_aa == val[1]:
      result = key
      break
  return result
    
def compute_hydrophobicity_character(aa_sequence, aa_symbol_size=3, span_size=23):
  aa_H = aa_hydrophobicity.aa_hydrophobicity() #get hydrophobicity of aminoacids
  protein_hydrophobicity = [] # the accumulator of result
  if aa_symbol_size == 3:
    aa_sequence_list = aa_sequence.split(sep="-") # aa seq in 3-letter with dash
  else:
    aa_sequence_list_1 = list(aa_sequence)
    aa_sequence_list = [convert_to_three_letter_aa(one_letter_aa) for\
      one_letter_aa in aa_sequence_list_1]
  if len(aa_sequence_list) >= span_size:
    for i in range(span_size//2, len(aa_sequence_list)-span_size//2, 1):
      temp = np.array( [aa_H[aa] for aa in aa_sequence_list[i-span_size//2:i+span_size//2] ] ) 
      protein_hydrophobicity.append( np.sum(temp) )
  return protein_hydrophobicity
    

if __name__ == "__main__":
  # x = compute_hydrophobicity_character(aa_sequence="Met-Ile-Leu-Thr-Hist", span_size=3)
  # print(x)
  # convert_to_three_letter_aa()
  x = compute_hydrophobicity_character(aa_sequence="MILYH",aa_symbol_size=1, span_size=5)
  print(len(x))

  
