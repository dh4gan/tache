# Written 9/10/14 by dh4gan
# Code reads in  output "eigen" file from sph_tidal_tensor
# Applies a user-defined threshold to the particles' tensor eigenvalues
# to classify the region they inhabit
# Classification is plotted in colour-coded format

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

# Colours of different classes
# Black= Void
# Green = filament
# Blue = sheet
# Red = cluster

colourchoice = ['black','green','blue','red']

# Read threshold value from command line or as argument

eigchar = raw_input("Analyse Velocity Shear or Tidal Tensor? (enter V or T): ")
individualchoice = raw_input("Plot each file separately? (y/n)")
threshold = input("What is the threshold for classification?")

prefix = "eigenvalues"+eigchar+"*"

filenames = filefinder.find_sorted_local_input_fileset(prefix)  
  
individualplots = False

if("y" in individualchoice or "Y" in individualchoice): individualplots=True

clusters = []
filaments = []
sheets = []
voids = []

for filename in filenames:

    # Read in eigenvalue file

    print "Reading eigenvalue file ",filename

    npart = read_eigenvalue_file.find_number_entries(filename)
    x,y,z,eigenpart,eigenvalues = read_eigenvalue_file.read_file(filename,npart)

    # fortran does rows and columns differently from python - switch them here
    eigenvalues = eigenvalues.transpose()

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

    clusters.append(np.size(classification[icluster)))
    filaments.append(np.size(classification[ifilament]))
    sheets.append(np.size(classification[isheet]))
    voids.append(np.size(classification[ivoid]))

    print "File: ",filename    
    print "Clusters: ",clusters[-1]
    print "Filaments: ",filaments[-1]
    print "Sheets: ", sheets[-1]
    print "Voids: ", voids[-1]
    print "-------"

    
    # Now plot this data

    if(individualplots):        

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
        
        outputfile = "classify_"+filename+"_"+str(threshold)+".png"
        print "Plotting to file ",outputfile
         
        plt.savefig(outputfile,format="png")

# Now plot the number of each class versus time

cumufig = plt.figure()
ax = cumufig.add_subplot(111)

ax.plot(np.log10(clusters), label='Clusters')
ax.plot(np.log10(sheets), label='Sheets')
ax.plot(np.log10(filaments), label='Filaments')
ax.plot(np.log10(voids), label='Voids')
ax.legend()
plt.show()

outputfile = "classes_vs_time.png"
print "Plotting to file ",outputfile
         
cumufig.savefig(outputfile,format="png")


print "Plotting complete"

