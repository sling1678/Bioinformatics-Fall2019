import numpy as np

OUTER_SIZE = 10
INNER_SIZE = 5

"""
Trapezoid rule is a more sophisticated implementation of the sliding window technique for calculating the hydrophobisity
of an amino acid.  It is implemented by a series of weights and two window sizes, the outer window size and the inner
window size.  The weights in the outer window on each side of the inner window increment up to the weight of the inner
window, whereas the inner window is of constant weight.
This function computes the weights for the position-weighted hydrophobicity average of the trapezoid rule.
Arguments:
  outer_size : int, number of residues to each side of center residue for the total window (i.e. half of total window size - 1)
  inner_size : int, number of residues to each side of center residue for the more narrow, inner window (i.e. half of inner window size - 1)
Returns:
  weights : List of weights to be used in calculation of hydrophibicity

Ref: J. Mol. Bid. (1992) 225, 487-494 
Membrane Protein Structure Prediction
Hydrophobicity Analysis and the Positive-inside Rule
Gunnar von Heijne 
"""
def compute_location_weights():
    # norm is a normalizing value used to compute outer window weights
    norm = (1 + OUTER_SIZE) ** 2 - INNER_SIZE ** 2

    # Computes weights for the outer window portion which is to the left of the inner window
    weights = [i / norm for i in range(1, OUTER_SIZE - INNER_SIZE + 2)]

    # Builds weighted values for inner window onto weight list.  Inner window weights are constant
    inner_weight = (OUTER_SIZE - INNER_SIZE + 1) / norm
    for i in range(OUTER_SIZE - INNER_SIZE + 2, OUTER_SIZE + INNER_SIZE + 1):
        weights.append(inner_weight)

    # Builds outer window weights for the portion to the right of the inner window by reversing the order of the
    # outer window weights on the left of the inner window
    for i in range(1, OUTER_SIZE - INNER_SIZE + 2):
        weights.append(weights[OUTER_SIZE - INNER_SIZE + 1 - i])

    return weights


def von_heijne_scale():
    "Ref: Von Heijne, J. Mol. Biol, 1992"
    hydrophobicities = {
        "A": 0.267, "C": 1.806, "D": -2.303, "E": -2.442, "F": 0.427, "G": 0.160,
        "H": -2.189, "I": 0.971, "K": -2.996, "L": 0.623, "M": 0.136, "N": -1.988,
        "P": -0.451, "Q": -1.814, "R": -2.749, "S": -0.119, "T": -0.083,
        "V": 0.721, "W": -0.875, "Y": -0.386,
    }
    return hydrophobicities

def aa_hydrophobicity():
    hydrophobicity = {
        "I"	: (4.5, 0.31,	-0.60,	-1.56,	1.97),
        "V"	: (4.2,	-0.07,	-0.31,	-0.78,	1.46),
        "L"	: (3.8,	0.56,	-0.55,	-1.81,	1.82),
        "F"	: (2.8,	1.13,	-0.32,	-2.20,	1.98),
        "C"	: (2.5,	0.24,	-0.13,	0.49,	-0.30),
        "M"	: (1.9,	0.23,	-0.10,	-0.76,	1.40),
        "A"	: (1.8,	-0.17,	0.11,	0.0,	0.38),
        "G"	: (-0.4,	-0.01,	0.74,	1.72,	-0.19),
        "T"	: (-0.7,	-0.14,	0.52,	1.78,	-0.32),
        "S"	: (-0.8,	-0.13,	0.84,	1.83,	-0.53),
        "W"	: (-0.9,	1.85,	0.30,	-0.38,	1.53),
        "Y"	: (-1.3,	0.94,	0.68,	-1.09,	0.49),
        "P"	: (-1.6,	-0.45,	2.23,	-1.52,	-1.44),
        "H"	: (-3.2,	-0.96,	2.06,	4.76,	-1.44),
        "E"	: (-3.5,	-2.02,	2.68,	1.64,	-2.90),
        "Q"	: (-3.5,	-0.58,	2.36,	3.01,	-1.84),
        "D"	: (-3.5,	-1.23,	3.49,	2.95,	-3.27),
        "N"	: (-3.5,	-0.42,	2.05,	3.47,	-1.62),
        "K"	: (-3.9,	-0.99,	2.71,	5.39,	-3.46),
        "R"	: (-4.5,	-0.81,	2.58,	3.71,	-2.57),
    }
    return hydrophobicity


"""
Uses hydrophobicity scale to build hydrophobicity profile using weighted averages calculated by the trapezoid rule

Arguments:
  aa_sequence : amino acid sequence

Returns:
  hp : hydrophobicity profile as numpy array
"""
def build_hydrophobicity_profile(aa_sequence):
    #hydrophobicities = von_heijne_scale()
    #hp_vals = np.array([hydrophobicities[aa] for aa in aa_sequence])

    hydrophobicities = aa_hydrophobicity()
    hp_vals = np.array([hydrophobicities[aa][0] for aa in aa_sequence])

    hp = []

    weights = np.array(compute_location_weights())

    if len(aa_sequence) >= OUTER_SIZE:
        for i in range(0, len(hp_vals) - weights.size + 1):
            hp.append(np.sum(hp_vals[i:i + weights.size].flatten() * weights))

    return hp


"""
Analyze each string of hydrophibic values and identify regions which are 'certain' transmembrane regions and which are
'putative' (aka, may or may not be transmembrane regions).

In von Heijne's paper, the lower_cutoff=0.5 and upper_cutoff=1.0, and we are supposed to work with 10 units on each side of the max.

Arguments:
    hp : weighted hydrophobic values (size = protein_size minus widow_size )
    upper_cutoff : (default=1.0) for definite inside membrane
    lower_cutoff : (default=0.5) for putative inside membrane
Returns:
    result : list of 'x' meaning undecided, 'M' definitely inside membrane, 'P' putatively inside membrane
"""
def analyze_hydrophobicity_profile(hp, upper_cutoff=1, lower_cutoff=.5):
    # Resultant array
    result_string = ['x' for i in range(len(hp) + 2 * OUTER_SIZE)]

    # Give each hydrophobic value an index corresponding to their position in the amino acid chain
    hp_with_index = [[hp[i], i] for i in range(len(hp))]

    # Sort hydrophobicity chain by hydrophobic value, order highest to lowest
    hp_with_index_sorted = sorted(hp_with_index, reverse=True)

    # Remove any hydrophobic acid from the list with a value less than the lower cutoff
    hp_with_index_relevant = [x for x in hp_with_index_sorted if x[0] >= lower_cutoff]

    # Work through the list, starting with the highest hydrophobic value.  If the window surrounding the next highest
    # value doesn't intersect with the window of a value already in our list, add that value to our list.
    already_done_set = set()
    hp_with_index_relevant2 = []
    for item in hp_with_index_relevant:
        temp_set = set(range(item[1] - OUTER_SIZE, item[1] + OUTER_SIZE))
        if already_done_set.intersection(temp_set):
            pass
        else:
            already_done_set = already_done_set | temp_set
            hp_with_index_relevant2.append(item)

    # Recall that our hydrophobic profile is (OUTER_SIZE * 2) smaller than the original amino acid sequence.  This is
    # because when we built our hyrdrophobic profile, the amino acids within OUTER_SIZE of either edge were left out of
    # the profile becuase they couldn't be properly measured by our sliding window technique.
    # Thus, we now increment each index with the value OUTER_SIZE, to get the index of each result as it relates to the
    # original amino acid sequence.
    hp_with_true_indices = [[x[0], x[1] + OUTER_SIZE] for x in hp_with_index_relevant2]

    # Sort the certain and potential hydrophobic regions by index positions (i.e. order in which they appear in the
    # amino acid sequence
    hp_with_true_indices = sorted(hp_with_true_indices, key = lambda x: x[1])

    # Replace the hydrophobic value of the critical residue with a simple marker indicating certain (M) versus
    # putative (P)
    # Replace index of critical residue with index of first residue in window
    # Add a second index indicating the last residue in the window.
    for item in hp_with_true_indices:
        if item[0] >= upper_cutoff:
            item[0] = 'M'
        elif item[0] >= lower_cutoff:
            item[0] = 'P'
        item[1] -= OUTER_SIZE
        item.append(item[1] + 2 * OUTER_SIZE)


    # Modify the resultant string array to show putative and certain transmembrane regions
    for item in hp_with_true_indices:
        if item[0] == 'M':
            for i in range(item[1], item[2]):
                result_string[i] = "M"
        elif item[0] == 'P':
            for i in range(item[1], item[2]):
                result_string[i] = "P"

    # Convert resultant string array to string
    result_string = ''.join(result_string)

    return result_string, hp_with_true_indices


def print_predictions(aa_sequence, result, line_len=50):
    cur_idx = 0
    while True:
        if cur_idx + line_len < len(aa_sequence):
            print(result[cur_idx:cur_idx + line_len])
            print(aa_sequence[cur_idx:cur_idx + line_len])
            cur_idx += line_len
        else:
            print(result[cur_idx:])
            print(aa_sequence[cur_idx:])
            break


if __name__ == "__main__":
    # aa_sequence="F"*7+"I"*7+"C"*7+"A"*7 # test sequence
    aa_sequence = "C" * 21
    hp = build_hydrophobicity_profile(aa_sequence)

    result = analyze_hydrophobicity_profile(hp)
    print(result)

def analyze_sequence(aminoacid_seq):
    hp = build_hydrophobicity_profile(aminoacid_seq)
    result = analyze_hydrophobicity_profile(hp)
    """
    print(f"length of amino acid sequence {len(aminoacid_seq)}")
    print(f"amino acid sequence{aminoacid_seq}")
    print(f"length of hydrophobicity profile {len(hp)}")
    print(f"hydrophobicity profile {hp}")
    print(f"length of result {len(result)}")
    print(f"result {result}")
    """
    return result
