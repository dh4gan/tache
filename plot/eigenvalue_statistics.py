# Written 9/10/14 by dh4gan
# Code reads in output eigenvalue file from tache
# Computes statistics

import numpy as np
import matplotlib.pyplot as plt
import io_tache as io

# Read in inputs from command line

filename = ff.find_local_input_files('eigenvalues*')
threshold = input("What is the threshold for classification? ")
  
# Read in eigenvalue file
print "Reading eigenvalue file ", filename
npart,x,y,z,eigenpart,eigenvalues = io.read_eigenvalue_file(filename)

print np.amax(eigenvalues),np.amin(eigenvalues)

# Calculate the trace for each simulation element

trace = np.zeros(npart)

for i in range(npart):
    for j in range(3):
        trace[i] = trace[i]+ eigenvalues[i,j]

normedeigenvalues = eigenvalues.copy()

for i in range(npart):
    if(trace[i]>0.0):
        normedeigenvalues[i,:] = normedeigenvalues[i,:]/trace[i]
    else:
        normedeigenvalues[i,:] = 0.0

# Make a histogram of the eigenvalues

alleigenvalues = eigenvalues.flatten()

fig1 = plt.figure(1)
ax = fig1.add_subplot(111)
ax.hist(alleigenvalues, bins=100, normed=True, log=True)
plt.show()


# Make a histogram of the traces

fig1 = plt.figure(1)
ax = fig1.add_subplot(111)
ax.hist(trace, bins=100, normed=True, log=True)
plt.show()






