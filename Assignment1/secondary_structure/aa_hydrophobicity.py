def aa_hydrophobicity():
  """
  Key : Residue Type	
  Value : (kdHydrophobicity(a),	wwHydrophobicity(b), hhHydrophobicity(c),
  	mfHydrophobicity(d),	ttHydrophobicity(e))
  Ref: https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/hydrophob.html
  (a) A simple method for displaying the hydropathic character of a protein. Kyte J, Doolittle RF. J Mol Biol. 1982 May 5;157(1):105-32.

  (b) Experimentally determined hydrophobicity scale for proteins at membrane interfaces. Wimley WC, White SH. Nat Struct Biol. 1996 Oct;3(10):842-8. Attribute assignment file wwHydrophobicity.txt.

  (c) Recognition of transmembrane helices by the endoplasmic reticulum translocon. Hessa T, Kim H, Bihlmaier K, Lundin C, Boekel J, Andersson H, Nilsson I, White SH, von Heijne G. Nature. 2005 Jan 27;433(7024):377-81, supplementary data. Attribute assignment file hhHydrophobicity.txt. In this scale, negative values indicate greater hydrophobicity.

  (d) Side-chain hydrophobicity scale derived from transmembrane protein folding into lipid bilayers. Moon CP, Fleming KG. Proc Natl Acad Sci USA. 2011 Jun 21;108(25):10174-7, supplementary data. Attribute assignment file mfHydrophobicity.txt. In this scale, negative values indicate greater hydrophobicity.

  (e) An amino acid “transmembrane tendency” scale that approaches the theoretical limit to accuracy for prediction of transmembrane helices: relationship to biological hydrophobicity. Zhao G, London E. Protein Sci. 2006 Aug;15(8):1987-2001. Attribute assignment file ttHydrophobicity.txt (contributed by Shyam M. Saladi).

  """  
  hydrophobicity = {
    "Ile"	: (4.5,	0.31,	-0.60,	-1.56,	1.97),
    "Val"	: (4.2,	-0.07,	-0.31,	-0.78,	1.46),
    "Leu"	: (3.8,	0.56,	-0.55,	-1.81,	1.82),
    "Phe"	: (	2.8,	1.13,	-0.32,	-2.20,	1.98),
    "Cys"	: (	2.5,	0.24,	-0.13,	0.49,	-0.30),
    "Met"	: (	1.9,	0.23,	-0.10,	-0.76,	1.40),
    "Ala"	: (	1.8,	-0.17,	0.11,	0.0,	0.38),
    "Gly"	: (	-0.4,	-0.01,	0.74,	1.72,	-0.19),
    "Thr"	: (	-0.7,	-0.14,	0.52,	1.78,	-0.32),
    "Ser"	: (	-0.8,	-0.13,	0.84,	1.83,	-0.53),
    "Trp"	: (	-0.9,	1.85,	0.30,	-0.38,	1.53),
    "Tyr"	: (	-1.3,	0.94,	0.68,	-1.09,	0.49),
    "Pro"	: (	-1.6,	-0.45,	2.23,	-1.52,	-1.44),
    "His"	: (	-3.2,	-0.96,	2.06,	4.76,	-1.44),
    "Glu"	: (	-3.5,	-2.02,	2.68,	1.64,	-2.90),
    "Gln"	: (	-3.5,	-0.58,	2.36,	3.01,	-1.84),
    "Asp"	: (	-3.5,	-1.23,	3.49,	2.95,	-3.27),
    "Asn"	: (	-3.5,	-0.42,	2.05,	3.47,	-1.62),
    "Lys"	: (	-3.9,	-0.99,	2.71,	5.39,	-3.46),
    "Arg"	: (	-4.5,	-0.81,	2.58,	3.71,	-2.57)
    }
  return hydrophobicity


if __name__ == "__main__":
  H = aa_hydrophobicity()
  print(H["Ile"][0]) #OK
