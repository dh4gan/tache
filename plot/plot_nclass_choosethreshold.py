# Written 9/10/14 by dh4gan
# Code reads in output "eigenvalues" file from tache
# Applies a user-defined threshold to the elements' tensor eigenvalues
# to classify the region they inhabit
# Classification is plotted via hexbin

import numpy as np
import matplotlib.pyplot as plt
import io_tache as io
import filefinder as ff
import read_eigenvalue_file
from sys import argv


filename = ff.find_local_input_files('eigenvalues*')
threshold = input("What is the threshold for classification? ")
  
# Read in eigenvalue file
print "Reading eigenvalue file ", filename
npart,x,y,z,eigenpart,eigenvalues = io.read_eigenvalue_file(filename)



#npart = read_eigenvalue_file.find_number_entries(filename)
#x,y,z,eigenpart,eigenvalues = read_eigenvalue_file.read_file(filename,npart)


# fortran does rows and columns differently from python - switch them here
#eigenvalues = eigenvalues.transpose()

tryagain = 'y'

# Now enter a do while loop until user is happy with result

while(tryagain!='n'):

    print "Classifying eigenvalues according to threshold ", threshold    
 
    classification = io.classify_all_eigenvalues(eigenvalues,npart,threshold)

    # Give some statistics on the classification

    nclusters,nfilaments,nsheets,nvoids = io.print_class_counts(classification,filename)

    # Plot classifications binned in x,y
    fig1 = io.plot_all_classes_xy(x,y,filename,classification,threshold)

    plt.show()
    
    # Check if user wants to try again
    
    tryagain = raw_input("Do you want to try a different threshold (y/n)?")
    
    if(tryagain!='n'):
        threshold = input("Enter a new threshold: ")
        
# End While Loop

# Once user is happy with the threshold, output this data to file

outputfile = "classified_"+filename+"_"+str(threshold)+".png"
fig1.savefig(outputfile,format="png")
