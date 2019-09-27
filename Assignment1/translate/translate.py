from translate import genetic_code

"""
Translates gene sequence into amino acid sequence
Assumptions:
    given mRNA starts with the start codon
Parameters:
    mRNA : string consisting of characters: AUGC 
Returns:
    aminoacid_sequence : one letter codings of amino acid sequence formed by mRNA; places '*' where there is codon error.
    nucleotide_errors : List of tuples from the mRNA sequence where the 3-nucleotide tuple is not in the dictionary of 
        amino acid codings.
"""


def translate_simple(mRNA):
    # Store the dictionary mapping mRNA tuple to their corresponding amino acid codings
    code = genetic_code.genetic_code()

    # Translate mRNA sequence into amino acid sequence
    aminoacid_sequence = ""
    nucleotide_errors = []
    for i in range(0, len(mRNA), 3):
        try:
            aminoacid_sequence += code[mRNA[i:i + 3]]
            if code[mRNA[i:i + 3]] == '_':
                # Stop translation when stop codon is reached
                break
        except KeyError:
            # Handle amino acid sequence errors
            aminoacid_sequence += "*"
            nucleotide_errors.append((i + 1, mRNA[i:i + 3]))

    return aminoacid_sequence[:-1], nucleotide_errors
