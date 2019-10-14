from deletion_insertion_and_substitution_costs import deletion_insertion_and_substitution_costs

def min_edit_distance(X, Y):
    """
    Parameters:
        X : str, (source) string to be changed to Y
        Y : str, target string
    Returns:
        D : edit distance matrix
    """
    
    deletion_costs, insertion_costs,\
         substitution_costs = deletion_insertion_and_substitution_costs()
    # Initialize
    N = len(X)
    M = len(Y)
    D = [[0 for j in range(M+1)] for i in range(N+1)]
    path = [[None for j in range(M+1)] for i in range(N+1)]
    for i in range(1,N+1):
        xi = i-1
        D[i][0] = D[i-1][0] + deletion_costs[X[xi]]
    for j in range(1, M+1):
        yj = j-1
        D[0][j] = D[0][j-1] + insertion_costs[Y[yj]]

    # Recursion
    for i in range(1, N+1):
        xi = i-1
        del_key = X[xi]
        for j in range(1, M+1):
            yj = j-1
            ins_key = Y[yj]
            sub_key = X[xi]+"-by-"+Y[yj]
            D[i][j] = min(
                D[i-1][j] + deletion_costs[ del_key ],
                D[i][j-1] + insertion_costs[ ins_key ],
                D[i-1][j-1] + substitution_costs[ sub_key ] 
            )
            if D[i][j] == D[i-1][j] + deletion_costs[ del_key ]:
                path[i][j] = "del " + del_key + ', '
            elif D[i][j] == D[i][j-1] + insertion_costs[ ins_key ]:
                path[i][j] = "ins " + ins_key + ', '
            else:
                path[i][j] = "sub " + sub_key + ', '

    return D, path

def print_path(path):
    i = len(path)-1
    j = len(path[0])-1
    last_to_first_path = [path[i][j]]
    while i > 1 and j > 1:
        if path[i-1][j-1] <= path[i][j]:
            i = i - 1
            j = j - 1
            last_to_first_path.append(path[i-1][j-1])
        elif path[i-1][j] <= path[i][j]:
            i = i - 1
            last_to_first_path.append(path[i-1][j])
        else:
            j = j - 1
            last_to_first_path.append(path[i][j-1])
    print(  last_to_first_path[::-1] )
    return  



if __name__ == "__main__":
    X = "ATGCA"
    Y = "ATGCA"
    D, path = min_edit_distance(X, Y)
    print(D)
    print( X, Y)
    print_path(path)