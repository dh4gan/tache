# Written 9/10/14 by dh4gan
# Code reads in multiple output "eigenvalues" file from tache
# Applies a user-defined threshold to the particles' tensor eigenvalues
# to classify the region they inhabit

import numpy as np
import matplotlib.pyplot as plt
import io_tache as io
import filefinder

# Check tensor file type, whether to make plots for every file, and threshold
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

    npart,x,y,z,eigenpart,eigenvalues = io.read_eigenvalue_file(filename)

    # Classify eigenvalues
    classification = io.classify_all_eigenvalues(eigenvalues,npart,threshold)

    # Give some statistics on the classification
    nclusters,nfilaments,nsheets,nvoids = io.print_class_counts(classification,filename)
    
    clusters.append(nclusters)
    filaments.append(nfilaments)
    sheets.append(nsheets)
    voids.append(nvoids)

    # Now plot this data if wished

    if(individualplots):      
        fig1 = io.plot_all_classes_xy(x,y,filename,classification,threshold)

        outputfile = "classify_"+filename+"_"+str(threshold)+".png"
        print "Plotting to file ",outputfile
         
        fig1.savefig(outputfile,format="png")
        plt.close(fig1)

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

