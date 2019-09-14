def genetic_code():
  """
  One and three-letter codes for aminoacids
  """
  one_letter_code = \
         {'UCA' : 'S',    # Serine
          'UCC' : 'S',    # Serine
          'UCG' : 'S',    # Serine
          'UCU' : 'S',    # Serine
          'UUC' : 'F',    # Phenylalanine
          'UUU' : 'F',    # Phenylalanine
          'UUA' : 'L',    # Leucine
          'UUG' : 'L',    # Leucine
          'UAC' : 'Y',    # Tyrosine
          'UAU' : 'Y',    # Tyrosine
          'UAA' : '_',    # Stop
          'UAG' : '_',    # Stop
          'UGC' : 'C',    # Cysteine
          'UGU' : 'C',    # Cysteine
          'UGA' : '_',    # Stop
          'UGG' : 'W',    # Tryptophan
          'CUA' : 'L',    # Leucine
          'CUC' : 'L',    # Leucine
          'CUG' : 'L',    # Leucine
          'CUU' : 'L',    # Leucine
          'CCA' : 'P',    # Proline
          'CCC' : 'P',    # Proline
          'CCG' : 'P',    # Proline
          'CCU' : 'P',    # Proline
          'CAC' : 'H',    # Histidine
          'CAU' : 'H',    # Histidine
          'CAA' : 'Q',    # Glutamine
          'CAG' : 'Q',    # Glutamine
          'CGA' : 'R',    # Arginine
          'CGC' : 'R',    # Arginine
          'CGG' : 'R',    # Arginine
          'CGU' : 'R',    # Arginine
          'AUA' : 'I',    # Isoleucine
          'AUC' : 'I',    # Isoleucine
          'AUU' : 'I',    # Isoleucine
          'AUG' : 'M',    # Methionine
          'ACA' : 'T',    # Threonine
          'ACC' : 'T',    # Threonine
          'ACG' : 'T',    # Threonine
          'ACU' : 'T',    # Threonine
          'AAC' : 'N',    # Asparagine
          'AAU' : 'N',    # Asparagine
          'AAA' : 'K',    # Lysine
          'AAG' : 'K',    # Lysine
          'AGC' : 'S',    # Serine
          'AGU' : 'S',    # Serine
          'AGA' : 'R',    # Arginine
          'AGG' : 'R',    # Arginine          
          'GUA' : 'V',    # Valine
          'GUC' : 'V',    # Valine
          'GUG' : 'V',    # Valine
          'GUU' : 'V',    # Valine
          'GCA' : 'A',    # Alanine
          'GCC' : 'A',    # Alanine
          'GCG' : 'A',    # Alanine
          'GCU' : 'A',    # Alanine
          'GAC' : 'D',    # Aspartic Acid
          'GAU' : 'D',    # Aspartic Acid
          'GAA' : 'E',    # Glutamic Acid
          'GAG' : 'E',    # Glutamic Acid
          'GGA' : 'G',    # Glycine
          'GGC' : 'G',    # Glycine
          'GGG' : 'G',    # Glycine
          'GGU' : 'G'}    # Glycine

  three_letter_code = \
        { 'UCA' : 'Ser',    # Serine
          'UCC' : 'Ser',    # Serine
          'UCG' : 'Ser',    # Serine
          'UCU' : 'Ser',    # Serine
          'UUC' : 'Phe',    # Phenylalanine
          'UUU' : 'Phe',    # Phenylalanine
          'UUA' : 'Leu',    # Leucine
          'UUG' : 'Leu',    # Leucine
          'UAC' : 'Tyr',    # Tyrosine
          'UAU' : 'Tyr',    # Tyrosine
          'UAA' : '___',    # Stop
          'UAG' : '___',    # Stop
          'UGC' : 'Cys',    # Cysteine
          'UGU' : 'Cys',    # Cysteine
          'UGA' : '___',    # Stop
          'UGG' : 'Trp',    # Tryptophan
          'CUA' : 'Leu',    # Leucine
          'CUC' : 'Leu',    # Leucine
          'CUG' : 'Leu',    # Leucine
          'CUU' : 'Leu',    # Leucine
          'CCA' : 'Pro',    # Proline
          'CCC' : 'Pro',    # Proline
          'CCG' : 'Pro',    # Proline
          'CCU' : 'Pro',    # Proline
          'CAC' : 'His',    # Histidine
          'CAU' : 'His',    # Histidine
          'CAA' : 'Gln',    # Glutamine
          'CAG' : 'Gln',    # Glutamine
          'CGA' : 'Arg',    # Arginine
          'CGC' : 'Arg',    # Arginine
          'CGG' : 'Arg',    # Arginine
          'CGU' : 'Arg',    # Arginine
          'AUA' : 'Ile',    # Isoleucine
          'AUC' : 'Ile',    # Isoleucine
          'AUU' : 'Ile',    # Isoleucine
          'AUG' : 'Met',    # Methionine
          'ACA' : 'Thr',    # Threonine
          'ACC' : 'Thr',    # Threonine
          'ACG' : 'Thr',    # Threonine
          'ACU' : 'Thr',    # Threonine
          'AAC' : 'Asn',    # Asparagine
          'AAU' : 'Asn',    # Asparagine
          'AAA' : 'Lys',    # Lysine
          'AAG' : 'Lys',    # Lysine
          'AGC' : 'Ser',    # Serine
          'AGU' : 'Ser',    # Serine
          'AGA' : 'Arg',    # Arginine
          'AGG' : 'Arg',    # Arginine          
          'GUA' : 'Val',    # Valine
          'GUC' : 'Val',    # Valine
          'GUG' : 'Val',    # Valine
          'GUU' : 'Val',    # Valine
          'GCA' : 'Ala',    # Alanine
          'GCC' : 'Ala',    # Alanine
          'GCG' : 'Ala',    # Alanine
          'GCU' : 'Ala',    # Alanine
          'GAC' : 'Asp',    # Aspartic Acid
          'GAU' : 'Asp',    # Aspartic Acid
          'GAA' : 'Glu',    # Glutamic Acid
          'GAG' : 'Glu',    # Glutamic Acid
          'GGA' : 'Gly',    # Glycine
          'GGC' : 'Gly',    # Glycine
          'GGG' : 'Gly',    # Glycine
          'GGU' : 'Gly'}    # Glycine


  return one_letter_code, three_letter_code
