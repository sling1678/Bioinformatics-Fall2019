def aminoacid_symbols():
  aa = {
    #Amino Acids with Hydrophobic Side Chain – Aliphatic
   "Ala" : ("Alanine", "A", "molecular_structures/Ala.gif"),
    "Ile" : ("Isoleucine", "I", "molecular_structures/Ile.gif"),
    "Leu" : ("Leucine",  "L", "molecular_structures/Leu.gif"),
    "Met" : ("Methionine", "M", "molecular_structures/Met.gif"),
    "Val" : ("Valine",  "V", "molecular_structures/Val.gif"),
    #Amino Acids with Hydrophobic Side Chain – Aromatic
         
    "Phe" : ("Phenylalanine", "F", "molecular_structures/Phe.gif"),
    "Trp" : ("Tryptophan", "W", "molecular_structures/Trp.gif"),
    "Tyr" : ("Tyrosine", "Y", "molecular_structures/Tyr.gif"),
    # Amino Acids with Polar Neutral Side Chains
    "Asn" : ("Asparagine", "N", "molecular_structures/Asn.gif"),
    "Cys" : ("Cysteine", "C", "molecular_structures/Cys.gif"),
    "Gln" : ("Glutamine", "Q", "molecular_structures/Gln.gif"), 
    "Ser" : ("Serine", "S", "molecular_structures/Ser.gif"), 
    "Thr" : ("Threonine", "T", "molecular_structures/Thr.gif"),
    #Amino Acids with Electrically Charged Side Chains – Acidic
    "Asp" : ("Aspartic acid", "D", "molecular_structures/Asp.gif"),
    "Glu" : ("Glutamic acid", "E", "molecular_structures/Glu.gif"),
    # Amino Acids with Electrically Charged Side Chains – Basic
    "Arg" : ("Arginine", "R", "molecular_structures/Arg.gif"),
    "His" : ("Histidine", "H", "molecular_structures/His.gif"),
    "Lys" : ("Lysine", "K", "molecular_structures/Lys.gif"),
    # Unique Amino Acids
    "Gly" : ("Glycine", "G", "molecular_structures/Gly.gif"),
    "Pro" : ("Proline", "P", "molecular_structures/Pro.gif"),

  }
  return aa
    

if __name__ == "__main__":
  aa  = aminoacid_symbols()
  print(aa["Lys"])

  