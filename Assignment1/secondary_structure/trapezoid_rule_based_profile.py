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

"""
Returns mapping of amino acid character codes to their approximate hydrophobicity mappings based on the Kyte, Doolittle
study.
"""
def KD_scale():
    hydrophobicity = {
        "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8,
        "G"	: -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3, "P": -1.6,
        "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
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
    hydrophobicities = KD_scale()
    hp_vals = np.array([hydrophobicities[aa] for aa in aa_sequence])

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
