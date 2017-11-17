# Written 4/3/15 by dh4gan
# Code reads in output eigenvalue files from tache
# Given a list of element IDs,
# it plots their classification as a function of time
#

import numpy as np
import matplotlib.pyplot as plt
from sys import argv
import filefinder as ff
import read_classification_file


labels = ['','Cluster', 'Filament','Sheet','Void']

# Read in elements to track from file
# File format should be a list of numbers, one per line

elementfile = raw_input("Enter the filename containing the element list: ")
eigchar = raw_input("Which tensor type: tidal (T) or velocityshear (V)? ")
eigchar = eigchar.upper()
threshold = input("What is the threshold for classification? ")

if(eigchar=='T'):
    print 'Using tidal tensor classification data'
elif(eigchar=='V'):
    print 'Using velocity shear tensor classification data'
else:
    print 'WARNING: unsure which classification to use - bad input: ',eigchar
    print 'Using tidal tensor data as default'
    eigchar = 'T'


print 'Reading element IDs from file ', elementfile
elementlist = np.genfromtxt(elementfile)
nelements = len(elementlist)

print 'Particle list read: total element number -', nelements 
print elementlist

# Find all files

filelist = ff.find_sorted_local_input_fileset('eigenvalues'+eigchar+'+*')

nfiles = len(filelist)

print 'Reading from ',nfiles, ' files:'

time = np.zeros(nfiles)
elementclasses = np.zeros((nelements,nfiles))

ifile = 0

for filename in filelist:

    # Read in eigenvalue file
    print "Reading eigenvalue file ", filename
    npart,x,y,z,eigenpart,eigenvalues = io.read_eigenvalue_file(filename)

    # Classify
    print "Classifying eigenvalues according to threshold ", threshold    
    classification = io.classify_all_eigenvalues(eigenvalues,npart,threshold)
            

    # Record element data in tracking arrays
    for ielement in range(nelements):
        elementclasses[ielement,ifile] = classification[elementlist[ielement]]
        
    
    ifile= ifile+1
    
# Save element data to output files

for ielement in range(nelements):
    
    filename = 'tracked_classifications.'+str(int(elementlist[ielement]))
    filedata = elementclasses[ielement,:]
    np.savetxt(filename,filedata)
    

# Plot all element classifications against time

fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.set_yticks((0,1,2,3,4))
ax.set_ylim(0,4)
ax.set_yticklabels(labels)
ax.set_xlabel('Time')
ax.set_ylabel('Classification')

for ielement in range(nelements):
    ax.plot(time,elementclasses[ielement,:])
plt.show()

fig1.savefig('elementclasses_vs_time.png',format='png')


