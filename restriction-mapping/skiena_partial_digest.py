#####Skiena's practical solution of assembling DNA from restriction pieces#####
# DNA is split into pieces and size of the pieces are determined by gel electrophoresis
# If restriction enzymes cleave DNA at x1, x2, x3, ..., xn, we will get:
# sizes of x2-x1, x3-x1, x4-x1, .., xn-x1, x3-x2, x4-x2, ..., xn-x_(n-1).
# Assume x1=0, and x2, x3, ..., xn>0 as points on the x-axis with scale given by 
# number of nucleotides. Then, Skiena's idea was to start by x=[x0=0, xn] and then
# decide on other coordinates by deciding on coordinates that the next largest piece 
# would generate. The algorithm goes by the name, partial_digest.
################################################################################

# Helper functions
import collections # for testing identity of two multisets
def is_subset(A, B):
    for a in A:
        if a not in B:
            return False
    return True
def remove_subset_elements(A, B):
    for a in A:
        B.remove(a)
    return B
def add_subset_elements(A, B):
    B.extend(A)
    return B
def is_equal_multi_set(A,B):
    return collections.Counter(A)==collections.Counter(B)

# The main function that decides one step of algorithm
def update_fragments_and_cleavage_sites(fragments, cleavage_sites, dna_length):
    """
    Removes fragments based on which assumption of coordinate of the cleavage site is the correct one.
    Arguments:
        fragments : list of fragment sizes remaining at this stage
        cleavage_sites: list of coordinates of the currently recognized cleavage sites
        dna_length : integer, the largest fragment in the original list of fragments
    Returns:
        updated fragments and cleavage_sites
    """
    frags = fragments.copy() # make a local shallow copy to work with
    cs = cleavage_sites.copy() # make a local shallow copy to work with

    frag_size_to_decide = max(frags)
    x_from_origin = frag_size_to_decide # if cleaved with x=0 and this x
    x_from_the_end = dna_length - x_from_origin # if cleaved from this x and x=xn

    # the multisets of fragments that would be generated in the case of the two 
    # possibilities
    deltaX = [abs(x-frag_size_to_decide) for x in cs ]
    deltaX_alt = [ abs(x_from_the_end - x) for x in cleavage_sites ]

    # Choose smallest x if the two cleavages give us the same fragments - just a choice
    # for the sake of consistency
    if is_equal_multi_set(deltaX, deltaX_alt):
        x_from_origin = min(x_from_the_end, x_from_origin) # pick the smaller

    # The logic of updating fragments and cleavage_sites
    if is_subset(deltaX, frags):
        cs.append(frag_size_to_decide)
        for x in deltaX:
            frags.remove(x)
    else:        
        if is_subset(deltaX_alt, frags):
            cs.append(x_from_the_end)
            for x in deltaX_alt:
                frags.remove(x)
    return frags, cs        

def partial_digest(fragments):
    """
    The main function. Given fragments, finds the cleavage sites.
    """
    remaining_frags = fragments.copy() # the local copy to be used to see if there are changes
    dna_length = max(fragments)
    fragments.remove(dna_length)
    cleavage_sites=[0, dna_length] # initialized with the two ends
    
    while fragments: 
        if sorted(remaining_frags) == sorted(fragments): # remaining fragments are inconsistent
            break
        fragments, cleavage_sites = \
            update_fragments_and_cleavage_sites(fragments, cleavage_sites, dna_length)
    return sorted( cleavage_sites ), sorted(fragments)

def main(fragments):
    cleavage_sites, fragments = partial_digest(fragments)
    if len(fragments) != 0:
        print( "Following fragments appear to come from some other dna, ", fragments )
    return cleavage_sites, fragments

def test_is_subset():
    B = [ 1, 2, 3, 4, 5]
    A = [2, 4]
    print(is_subset(A, B)) #OK

def test_remove_subset_elements():
    B = [ 1, 2, 3, 4, 5]
    A = [2, 4]    
    print(B, A, remove_subset_elements(A,B)) #OK

def test_add_subset_elements():
    B = [ 1, 2, 3, 4, 5]
    A = [2, 4]    
    print(B, A, add_subset_elements(A,B))  #OK  


def test_update_fragments_and_cleavage_sites():
    dna_length = 10
    fragments = [2, 3, 4, 6]
    cleavage_sites = [0, 10, 2, 7]
    print("fragments, cleavage_sites: ", fragments, cleavage_sites)
    fragments, cleavage_sites = update_fragments_and_cleavage_sites(fragments, cleavage_sites, dna_length)
    print("fragments, cleavage_sites: ", fragments, cleavage_sites)

def test_partial_digest():
    fragments = [2, 2, 3, 3, 4, 5, 6, 7, 8, 10]
    cleavage_sites, remaining_fragments = partial_digest(fragments)
    print("cleavage_sites: ", cleavage_sites)
    print("remaining_fragments: ", remaining_fragments)

def test_main():

    fragments = [2, 2, 3, 3, 4, 5, 6, 7, 8, 10]
    cleavage_sites, remaining_fragments =main(fragments)
    print(cleavage_sites)
if __name__ == "__main__":   
    # cleavage_sites = partial_digest(fragments)
    # print(cleavage_sites)
    # test_is_subset()
    # test_remove_subset_elements()
    # test_add_subset_elements()
    # test_update_fragments_and_cleavage_sites()
    # test_partial_digest()
    test_main()




