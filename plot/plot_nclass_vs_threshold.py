# Written 9/10/14 by dh4gan
# Code reads in  output "eig" file from tache
# Applies multiple thresholds to the particles' tensor eigenvalues
# to classify the region they inhabit
# The code then outputs statistics for each threshold value

import numpy as np
import matplotlib.pyplot as plt
import classify_eigenvalues as classify
import read_eigenvalue_file
import filefinder
from sys import argv

# Column indices for different properties
xcol = 0
ycol = 1
zcol = 2
eigcols = range(3,5)

# Selection of thresholds

thresholdlist = np.logspace(-3, 3, num=100)

# Colours of different classes
# Black= Void
# Green = filament
# Blue = sheet
# Red = cluster

colourchoice = ['black','green','blue','red']

# Read threshold value from command line or as argument

if(len(argv)==1):
    filename = raw_input("What is the eigenvalue filename? ")
elif(len(argv)==2):
    filename = argv[1]
      
print "Finding Particle Number"
npart = read_eigenvalue_file.find_number_entries(filename)

print "Reading eigenvalue file ", filename
x,y,z,eigenpart,eigenvalues = read_eigenvalue_file.read_file(filename,npart)

# fortran does rows and columns differently from python - switch them here
eigenvalues = eigenvalues.transpose()

clusters = []
filaments = []
sheets = []
voids = []

thresholdcounter = 0
for threshold in thresholdlist:

    thresholdcounter = thresholdcounter +1
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

    clusters.append(np.size(classification[classification[:]==0]))
    filaments.append(np.size(classification[classification[:]==1]))
    sheets.append(np.size(classification[classification[:]==2]))
    voids.append(np.size(classification[classification[:]==3]))

    print "Threshold ",thresholdcounter, ": ", threshold
    print "Clusters: ",clusters[-1]
    print "Filaments: ",filaments[-1]
    print "Sheets: ", sheets[-1]
    print "Voids: ", voids[-1]
    print "-------"

    
# Now plot this data

clusters[:] = [c/float(npart) for c in clusters]
filaments[:] = [c/float(npart) for c in filaments]
sheets[:] = [c/float(npart) for c in sheets]
voids[:] = [c/float(npart) for c in voids]

fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(1.0e-6,1.0e0)


if(np.sum(clusters)>0.0): ax.plot(thresholdlist,clusters, color='blue', label='Clusters')
if(np.sum(filaments)>0.0): ax.plot(thresholdlist,filaments,color='green', label='Filaments')
if(np.sum(sheets)>0.0): ax.plot(thresholdlist,sheets, color='red',label='Sheets')
if(np.sum(voids)>0.0): ax.plot(thresholdlist,voids, color='cyan',label='Voids')
ax.set_xlabel('Threshold eigenvalue (dimensionless)')
ax.set_ylabel('Fraction of Simulation')
ax.legend(loc = 'lower right')
plt.show()

fig1.savefig('multithreshold_'+filename+'.png', format='png')
    
print "Plotting complete"
