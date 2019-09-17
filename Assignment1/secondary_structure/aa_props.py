import pandas as pd 
import numpy as np
def aa_props(aa_prop_csv_data_file):
  """
  Reference:
  https://www.sigmaaldrich.com/life-science/metabolomics/learning-center/amino-acid-reference-chart.html
  1 pKa is the negative of the logarithm of the dissociation constant for the -COOH group.
  2 pKb is the negative of the logarithm of the dissociation constant for the -NH3 group.
  3 pKx is the negative of the logarithm of the dissociation constant for any other group in the molecule.
  4 pl is the pH at the isoelectric point.
  Reference: D.R. Lide, Handbook of Chemistry and Physics, 72nd Edition, CRC Press, Boca Raton, FL, 1991.

  The table has following props:
  Name	3-Letter Symbol	1-Letter Symbol	Molecular Weight	Molecular
   Formula	Residue Formula	Residue Weight (-H2O)	pKa1	pKb2	pKx3	pI4

  """
  aa_dataframe = pd.read_csv(aa_prop_csv_data_file)
  return aa_dataframe

if __name__ == "__main__":
  df = aa_props(aa_prop_csv_data_file="./aa_props_table.csv")
  print(df.columns)
  print(df[['Name', 'Charge', 'Hydrophobicity']])