from translate import genetic_code


"""
Translates gene sequence into amino acid sequence
Assumptions:
  given mRNA starts with a start codon.
Parameters:
  mRNA : string consisting of AUGC and multiple of 3
Returns:
  aminoacid_sequence : get one-letter sequence; places '*' where there is codon error.
  nucleotide_errors : List of tuple (index of mRNA sequence where a 3-nucleotide is not in the dictionary of codes
"""
def translate_simple(mRNA):
  # Codon to amino acid dictionary
  code = genetic_code.genetic_code()
  
  # Translate
  aminoacid_sequence = ""
  nucleotide_errors = []
  for i in range(0, len(mRNA), 3):
    try:
      aminoacid_sequence += code[mRNA[i:i+3]]
      if code[mRNA[i:i+3]] == '_':
        # Stop translation when stop codon is reached
        break
    except KeyError:
      # Handle amino acid sequence errors
      aminoacid_sequence += "*"
      nucleotide_errors.append((i+1, mRNA[i:i+3]))

  return aminoacid_sequence[:-1], nucleotide_errors
