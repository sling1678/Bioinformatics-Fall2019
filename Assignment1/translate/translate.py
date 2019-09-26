from translate import genetic_code


def translate_simple(mRNA):
  """
  Assumes given mRNA starts with a start codon.
  Parameters:
    mRNA : string consisting of AUGC and multiple of 3
    one_letter : set to True if want one-letter, default False
  Returns:
    aminoacid_sequence : get one-letter sequence; places '*' where there is codon error.
    nucleotide_errors : List of tuple (index of mRNA sequence where a 3-nucleotide is not in the dictionary of codes
  """
  code = genetic_code.genetic_code()
  
  assert(len(mRNA)%3 == 0), "mRNA must have multiples of 3 nucleotides"
  aminoacid_sequence = ""
  nucleotide_errors = []
  for i in range(0,len(mRNA), 3):
    try:
      aminoacid_sequence += code[mRNA[i:i+3]]
    except KeyError:
      aminoacid_sequence += "*"
      nucleotide_errors.append((i+1, mRNA[i:i+3]))
    if code[mRNA[i:i+3]] == '_':
      break

  return aminoacid_sequence[:-1], nucleotide_errors
