# Written 9/10/14 by dh4gan
# Code reads in output "eig" file from tache
# Applies a user-defined threshold to the particles' tensor eigenvalues
# to classify the region they inhabit
# Classification is plotted in colour-coded format

import numpy as np
import matplotlib.pyplot as plt
import classify_eigenvalues as classify
import read_eigenvalue_file
from sys import argv

# Column indices for different properties
xcol = 0
ycol = 1
zcol = 2
eigcols = range(3,6)

# Colours of different classes
# Black= Void
# Green = filament
# Blue = sheet
# Red = cluster

colourchoice = [[0,0,0],[0,1,0],[0,0,1],[1,0,0]]

# Read threshold value from command line or as argument

if(len(argv)==1):
    threshold = input("What is the threshold for classification? ")
    filename = raw_input("What is the eigenvalue filename? ")
elif(len(argv)==2):
    threshold = float(argv[1])
    filename = raw_input("What is the eigenvalue filename?")
elif(len(argv)>2):
    threshold = float(argv[1])
    filename = argv[2]
  
# Read in eigenvalue file

print "Finding Particle Number"

npart = read_eigenvalue_file.find_number_entries(filename)

print "Reading eigenvalue file ", filename

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
    colourpart = []

    for i in range(npart):
        eig = eigenvalues[i,:]
    
        classification[i] = classify.classify_eigenvalue(eig, threshold)
        colourpart.append(colourchoice[classification[i]])


    # Give some statistics on the classification

    colourpart = np.array(colourpart)

    clusters = classification[classification[:]==0]
    filaments = classification[classification[:]==1]
    sheets = classification[classification[:]==2]
    voids = classification[classification[:]==3]

    print "Clusters: ",np.size(clusters)
    print "Filaments: ",np.size(filaments)
    print "Sheets: ", np.size(sheets)
    print "Voids: ", np.size(voids)

    # Now plot this data

    xmin = np.amin(x)
    xmax = np.amax(x)

    ymin = np.amin(y)
    ymax = np.amax(y)

    #xmin = -300.0
    #xmax = 300.0

    #ymin = -300.0
    #ymax = 300.0

    markersize = 10.0

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
    axcl.scatter(x[classification[:]==0],y[classification[:]==0], c=colourpart[classification[:]==0], s=markersize)
    axcl.set_xlim(xmin,xmax)
    axcl.set_ylim(ymin,ymax)

    axfil.set_title("Filaments")
    axfil.scatter(x[classification[:]==1],y[classification[:]==1], c=colourpart[classification[:]==1], s=markersize)
    axfil.set_xlim(xmin,xmax)
    axfil.set_ylim(ymin,ymax)

    axsh.set_title("Sheets")
    axsh.scatter(x[classification[:]==2],y[classification[:]==2], c=colourpart[classification[:]==2], s = markersize)
    axsh.set_xlim(xmin,xmax)
    axsh.set_ylim(ymin,ymax)

    axvd.set_title("Voids")
    axvd.scatter(x[classification[:]==3],y[classification[:]==3], c=colourpart[classification[:]==3], s = markersize)
    axvd.set_xlim(xmin,xmax)
    axvd.set_ylim(ymin,ymax)

    plt.show()
    
    # Check if user wants to try again
    
    tryagain = raw_input("Do you want to try a different threshold (y/n)?")
    
    if(tryagain!='n'):
        threshold = input("Enter a new threshold: ")
        
# End While Loop

# Once user is happy with the threshold, output this data to file

outputfile = "classify_"+filename+"_"+str(threshold)+".png"
fig1.savefig(outputfile,format="png")
