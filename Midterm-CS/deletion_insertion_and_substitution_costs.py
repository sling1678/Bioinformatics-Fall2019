def deletion_insertion_and_substitution_costs():
    deletion_costs = {
        "A": 1, "T": 1, "G": 1, "C": 1
    }
    insertion_costs = {
        "A": 1, "T": 1, "G": 1, "C": 1
    }
    substitution_costs = {
        # A, G: purines; T, C: pyrimidines
        # Substituting Purines
        "A-by-G" : 1, "G-by-A" : 1,
        "A-by-T" : 2, "A-by-C" : 2, "G-by-T" : 2, "G-by-C" : 2,
        # Substituting Pyrimidines
        "T-by-C" : 1, "C-by-T" : 1,
        "T-by-G" : 2, "T-by-A" : 2, "C-by-G" : 2, "C-by-A" : 2,
        # Self subs
        "A-by-A" : 0, "G-by-G" : 0, "T-by-T" : 0, "C-by-C" : 0
    }
    return deletion_costs, insertion_costs, substitution_costs    