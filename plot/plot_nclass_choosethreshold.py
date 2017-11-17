# Written 9/10/14 by dh4gan
# Code reads in output "eigenvalues" file from tache
# Applies a user-defined threshold to the elements' tensor eigenvalues
# to classify the region they inhabit
# Classification is plotted via hexbin

import numpy as np
import matplotlib.pyplot as plt
import classify_eigenvalues as classify
import filefinder as ff
import read_eigenvalue_file
from sys import argv


filename = ff.find_local_input_files('eigenvalues*')
threshold = input("What is the threshold for classification? ")
  
# Read in eigenvalue file
print "Reading eigenvalue file ", filename

npart = read_eigenvalue_file.find_number_entries(filename)
x,y,z,eigenpart,eigenvalues = read_eigenvalue_file.read_file(filename,npart)


# fortran does rows and columns differently from python - switch them here
eigenvalues = eigenvalues.transpose()

tryagain = 'y'

# Now enter a do while loop until user is happy with result

while(tryagain!='n'):

    print "Classifying eigenvalues according to threshold ", threshold    
    # Classify eigenvalues
    # Function returns an integer iclass
    #     iclass = 0 --> cluster
    #     iclass = 1 --> filament
    #     iclass = 2 --> sheet
    #     iclass = 3 --> void

    classification = np.empty(npart, dtype="int")


    for i in range(npart):
        eig = eigenvalues[i,:]
    
        classification[i] = classify.classify_eigenvalue(eig, threshold)

    # Give some statistics on the classification

    icluster = classification[:]==0
    ifilament = classification[:]==1
    isheet = classification[:]==2
    ivoid = classification[:]==3

    clusters = classification[icluster]
    filaments = classification[ifilament]
    sheets = classification[isheet]
    voids = classification[ivoid]

    print "Clusters: ",np.size(clusters)
    print "Filaments: ",np.size(filaments)
    print "Sheets: ", np.size(sheets)
    print "Voids: ", np.size(voids)

    # Now plot this data

    xmin = np.amin(x)
    xmax = np.amax(x)

    ymin = np.amin(y)
    ymax = np.amax(y)

    fig1, ((axcl,axfil),(axsh,axvd)) = plt.subplots(2,2, sharex='col',sharey='row')

    # Decide on plot title using filename:
    # Files containing "T" have been calculated using the tidal tensor
    # Files containing "V" have been calculated using the velocity shear tensor
    
    plot_title = "Tensor Classification, Threshold ="+str(threshold)
    
    if "T" in filename:
        plot_title = "Tidal Tensor Classification, Threshold ="+str(threshold)
    
    if "V" in filename:
        plot_title = "Velocity Shear Tensor Classification, Threshold ="+str(threshold)

    plt.suptitle(plot_title)
    axcl.set_title("Clusters")
    axcl.hexbin(x[icluster],y[icluster],gridsize=500)
    axcl.set_xlim(xmin,xmax)
    axcl.set_ylim(ymin,ymax)

    axfil.set_title("Filaments")
    axfil.hexbin(x[ifilament],y[ifilament],gridsize=500)
    axfil.set_xlim(xmin,xmax)
    axfil.set_ylim(ymin,ymax)

    axsh.set_title("Sheets")
    axsh.hexbin(x[isheet],y[isheet],gridsize=500)
    axsh.set_xlim(xmin,xmax)
    axsh.set_ylim(ymin,ymax)

    axvd.set_title("Voids")
    axvd.hexbin(x[ivoid],y[ivoid],gridsize=500)
    axvd.set_xlim(xmin,xmax)
    axvd.set_ylim(ymin,ymax)

    plt.show()
    
    # Check if user wants to try again
    
    tryagain = raw_input("Do you want to try a different threshold (y/n)?")
    
    if(tryagain!='n'):
        threshold = input("Enter a new threshold: ")
        
# End While Loop

# Once user is happy with the threshold, output this data to file

outputfile = "classified_"+filename+"_"+str(threshold)+".png"
fig1.savefig(outputfile,format="png")
