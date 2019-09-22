from translate import genetic_code


def translate_simple(mRNA, one_letter=False):
  """
  Assumes given mRNA starts with a start codon.
  Parameters:
    mRNA : string consisting of AUGC and multiple of 3
    one_letter : set to True if want one-letter, default False
  Returns:
    aminoacid_sequence : default 3-letter symbols; with one-letter set to True, get one-letter sequence; places '*' where there is codon error.
    nucleotide_errors : List of tuple (index of mRNA sequence where a 3-nucleotide is not in the dictionary of codes
  """
  if one_letter:
    code, _ = genetic_code.genetic_code()
  else:
    _, code = genetic_code.genetic_code()
  
  assert(len(mRNA)%3 == 0), "mRNA must have multiples of 3 nucleotides"
  aminoacid_sequence = ""
  nucleotide_errors = []
  for i in range(0,len(mRNA), 3):
    try:
      aminoacid_sequence += code[mRNA[i:i+3]]
    except KeyError:
      aminoacid_sequence += "*"
      nucleotide_errors.append((i+1, mRNA[i:i+3]))
    if not one_letter:
      if i < len(mRNA)-3:
        aminoacid_sequence += "-"
  return aminoacid_sequence, nucleotide_errors
