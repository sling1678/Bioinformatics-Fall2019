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
    """norm is a normalizing value used to compute outer window weights"""
    norm = (1 + OUTER_SIZE) ** 2 - INNER_SIZE ** 2

    """Computes weights for the outer window portion which is to the left of the inner window"""
    weights = [i / norm for i in range(1, OUTER_SIZE - INNER_SIZE + 2)]

    """Builds weighted values for inner window onto weight list.  Inner window weights are constant"""
    inner_weight = (OUTER_SIZE - INNER_SIZE + 1) / norm
    for i in range(OUTER_SIZE - INNER_SIZE + 2, OUTER_SIZE + INNER_SIZE + 1):
        weights.append(inner_weight)

    """Builds outer window weights for the portion to the right of the inner window by reversing the order of the 
  outer window weights on the left of the inner window"""
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


"""
Uses hydrophobicity scale to build hydrophobicity profile using weighted averages calculated by the trapezoid rule

Arguments:
  aa_sequence : amino acid sequence

Returns:
  hp : hydrophobicity profile as numpy array

"""


def build_hydrophobicity_profile(aa_sequence):
    hydrophobicities = von_heijne_scale()
    hp_vals = np.array([hydrophobicities[aa] for aa in aa_sequence])

    hp = []

    weights = np.array(compute_location_weights())

    if len(aa_sequence) >= OUTER_SIZE:
        for i in range(0, len(hp_vals) - weights.size - 1):
            hp.append(np.sum(hp_vals[i:i + weights.size].flatten() * weights))
    return hp


"""
Classify the membrane segments as : types 'certain' and 'piutative'.
Procedure: Find highest hp value, then classify the 21 residues around it as certain or putative depending on the value at the peak, then look for next pean and repeat till less than 20 residues left or highest peak is less than lower_cutoff.

In von Heijne's paper, the lower_cutoff=0.5 and upper_cutoff=1.0, and we are supposed to work with 10 units on each side of the max.

Procedure:
    1. Get values and indices of all the local maxima points whose peak values are at least the lower_cutoff
    3. Order these (peak-value, protein-indice) using peak-value as key
    4. Remove lower peak values if they fall within plus-minus of the outer_size
    5. convert these indices to indices in the protein sequence
    6. Now process these values as follows:
        initialize result as 'x' for every aa in sequence
        if peak_value >= upper_cutoff, then label corresponding index-outer:index+outer+1 aminoacids as definite_membrane (M) - the 'x' will become 'M'
        elif peak_value >= lower_cutoff, then label corresponding index-outer:index+outer+1 aminoacids as putative_membrane (?) - the 'x' will become '?'
    7. return result
    8. Later print one line of prediction over one line of one-letter aa-sequence
    9. More advanced - process topology by positive-inside rule
Arguments:
    hp : weighted hydrophobic values (size = protein_size minus widow_size )
    outer_size : window size
    upper_cutoff : (default=1.0) for definite inside membrane
    lower_cutoff : (default=0.5) for putative inside membrane
    critical_num_residues : (default=21) for processing topologies - not used so far
Returns:
    result : list of 'x' meaning undecided, 'M' definitely inside membrane, 'P' putatively inside membrane
"""


def analyze_hydrophobicity_profile(hp, upper_cutoff=1.0, lower_cutoff=0.5, critical_num_residues=21):
    result = ['x' for i in range(len(hp) + 2 * OUTER_SIZE)]

    hp_with_index = [[hp[i], i] for i in range(len(hp))]
    hp_with_index_sorted = sorted(hp_with_index, reverse=True)

    hp_with_index_relevant = [x for x in hp_with_index_sorted if x[0] >= lower_cutoff]

    already_done_set = set()
    hp_with_index_relevant2 = []
    for item in hp_with_index_relevant:
        temp_set = set(range(item[1] - OUTER_SIZE, item[1] + OUTER_SIZE + 1, 1))
        if already_done_set.intersection(temp_set):
            pass
        else:
            already_done_set = already_done_set | temp_set
            hp_with_index_relevant2.append(item)

    hp_with_true_indices = [[x[0], x[1] + OUTER_SIZE] for x in hp_with_index_relevant2]

    for item in hp_with_true_indices:
        if item[0] >= upper_cutoff:
            for i in range(item[1] - OUTER_SIZE, item[1] + OUTER_SIZE + 1, 1):
                result[i] = "M"
        elif item[0] >= lower_cutoff:
            for i in range(item[1] - OUTER_SIZE, item[1] + OUTER_SIZE + 1, 1):
                result[i] = "P"

    result = ''.join(result)

    return result


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
    aa_sequence = "I" * 7 + "R" * 1 + "C" * 5 + "G" * 10 + "F" * 20 + "I" * 25  # test sequence
    hp = build_hydrophobicity_profile(aa_sequence)

    result = analyze_hydrophobicity_profile(hp)
    print(aa_sequence)
    print(len(aa_sequence))
    print(result)
    print(len(result))
    # print the output
    # print(len(result) == len(aa_sequence))
    # print_predictions(aa_sequence, result)
