# Written 9/10/14 by dh4gan
# Some useful functions for classifying eigenvalues and defining structure

def classify_eigenvalue(eigenvalues, threshold):
    '''Given 3 eigenvalues, and some threshold, returns an integer
    'iclass' corresponding to the number of eigenvalues below the threshold
    
    iclass = 0 --> clusters (3 +ve eigenvalues, 0 -ve)
    iclass = 1 --> filaments (2 +ve eigenvalues, 1 -ve)
    iclass = 2 --> sheet (1 +ve eigenvalues, 2 -ve)
    iclass = 3 --> voids (0 +ve eigenvalues, 3 -ve)
    
    '''
    
    iclass = 0
    
    for i in range(3):

        if(eigenvalues[i]<threshold):
            iclass +=1
            
    return int(iclass)
