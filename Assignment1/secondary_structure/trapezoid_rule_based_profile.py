import numpy as np 
import pandas as pd 
import os

def compute_location_weights(outer_size=10, inner_size=5):
  """
  Computes weights for position-weighted hydrophobicity average for trapezoid rule.
  Arguments:
    outer_size : int, the number of residues 2*outer_size+1 in wider window
    inner_size : int, the number of residues 2*inner_size+1 in narrow window
  Returns:
    weights : List of weigths to be used by residues in window of size outer_size

  Ref: J. Mol. Bid. (1992) 225, 487-494 
  Membrane Protein Structure Prediction
  Hydrophobicity Analysis and the Positive-inside Rule
  Gunnar von Heijne 
  """

  n= outer_size
  q=inner_size
  norm = (1 + n)**2 - q**2
  weights = [i/norm for i in range(1,n-q+2,1)]
  for i in range(n-q+2, n+q+1, 1):
    weights.append( (n-q+1)/norm )
  for i in range(n+q+1, 2*n+2,1):
    weights.append( (2*n+2-i)/norm )
  return weights

def get_aa_props(path_to_aa_props_table):
  aa_dataframe = pd.read_csv(path_to_aa_props_table)
  aa_dataframe = aa_dataframe.set_index("Symbol_3", drop = False)
  return aa_dataframe

def von_heijne_scale():
  "Ref: Von Heijne, J. Mol. Biol, 1992"
  hydrophobicities ={
    "A": 0.267, "C": 1.806, "D": -2.303, "E": -2.442, "F": 0.427, "G": 0.160,
    "H": -2.189, "I": 0.971, "K": -2.996, "L": 0.623, "M": 0.136, "N":-1.988,
    "P": -0.451, "Q": -1.814, "R": -2.749, "S": -0.119, "T": -0.083, 
    "V": 0.721, "W": -0.875, "Y": -0.386,
  }
  return hydrophobicities

def build_hydrophobicity_profile(aa_sequence_1_let, hydrophobicities, outer_size=3, inner_size=1):
  """
  Uses hydrophobicity scale to build hydrophobicity profile by weighted average over a trapezoid rule

  Arguments:
    aa_sequence_1_let : one-letter symbols for aminoacid sequence
    hydrophobicities : dictionary that maps one-letter symbols of aa to values
    outer_size : int, the number of residues 2*outer_size+1 in wider window
    inner_size : int, the number of residues 2*inner_size+1 in narrow window 
  Returns:
    hp : hydrophobicity profile as numpy array
  """
  hp_vals=np.array( [ hydrophobicities[aa] for aa in aa_sequence_1_let]  )
 
  hp = []
  weights = np.array( compute_location_weights(outer_size, inner_size) )
  if len(aa_sequence_1_let) < outer_size:
    pass
  else:
    hp = [ np.sum((hp_vals[i:i+2*outer_size+1]).flatten()*weights) for i in range(0, len(hp_vals)-2*outer_size,1) ]
  return hp, outer_size




def analyze_huydrophobicity_profile(hp, outer_size, lowest_cutoff=0.5, mid_cutoff=1.0, critical_num_residues=21):
  """
  Classify the membrane segments as : types 'certain' and 'piutative'.
  Procedure: Find highest hp value, then classify the 21 residues around it as certain or putative depending on the value at the peak, then look for next pean and repeat till less than 20 residues left or highest peak is less than lower_cutoff.

  In von Heijne's paper, the lower_cutoff=0.5 and upper_cutoff=1.0, and we are supposed to work with 10 units on each side of the max.

  Procedure:

    1. Get values and indices of all the local maxima points whose peak values are at lease the lower_cutoff
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
  """
  result = [] 
  # TODO

  




if __name__ == "__main__":
  dir_path = os.path.dirname(os.path.realpath(__file__))
  path_to_aa_props_table = dir_path+"/"+"aa_props_table.csv"
  aa_dataframe = get_aa_props(path_to_aa_props_table)

  aa_sequence="F"*7+"I"*7+"C"*7+"A"*7
  hydrophobicities = von_heijne_scale()
  hp, outer_size = build_hydrophobicity_profile(aa_sequence, hydrophobicities, outer_size=3, inner_size=1)
  print(hp)
  
  # analyze_huydrophobicity_profile(hp, outer_size, lowest_cutoff=-2.5, mid_cutoff=-1.7, critical_num_residues=21)