# Written 9/10/14 by dh4gan
# Code reads in  output "eigenvalues" file from tache
# Applies multiple thresholds to the particles' tensor eigenvalues
# to classify the region they inhabit
# The code then outputs statistics for each threshold value

import numpy as np
import matplotlib.pyplot as plt
import io_tache as io
import filefinder

# Selection of thresholds

thresholdlist = np.logspace(-3, 3, num=100)

filename = ff.find_local_input_files('eigenvalues*')
  
# Read in eigenvalue file
print "Reading eigenvalue file ", filename

npart,x,y,z,eigenpart,eigenvalues = io.read_eigenvalue_file.find_number_entries(filename)

clusters = []
filaments = []
sheets = []
voids = []

thresholdcounter = 0
for threshold in thresholdlist:

    thresholdcounter = thresholdcounter +1

    # Classify eigenvalues
    classification = io.classify_all_eigenvalues(eigenvalues,npart,threshold)

    # Give some statistics
    print "Threshold ",thresholdcounter, ": ", threshold
    nclusters,nfilaments,nsheets,nvoids = io.print_class_counts(classification,filename)

    clusters.append(nclusters)
    filaments.append(nfilaments)
    sheets.append(nsheets)
    voids.append(nvoids)
    
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


if(np.sum(clusters)>0.0): 
    ax.plot(thresholdlist,clusters, color='blue', label='Clusters')
if(np.sum(filaments)>0.0): 
    ax.plot(thresholdlist,filaments,color='green', label='Filaments')
if(np.sum(sheets)>0.0): 
    ax.plot(thresholdlist,sheets, color='red',label='Sheets')
if(np.sum(voids)>0.0): 
    ax.plot(thresholdlist,voids, color='cyan',label='Voids')

ax.set_xlabel('Threshold eigenvalue (dimensionless)')
ax.set_ylabel('Fraction of Simulation')
ax.legend(loc = 'lower right')
plt.show()

fig1.savefig('multithreshold_'+filename+'.png', format='png')
    
print "Plotting complete"
