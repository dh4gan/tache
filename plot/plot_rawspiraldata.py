# Written 02/02/17 by dh4gan
# Reads in the output from FORTRAN 90 code 'spiralfind'
# (List of spiral files, with x y z points)

# Plots the raw spiral data for quick inspection

# Runs an MCMC (Metropolis-Hastings) to find best fitting logarithmic spiral
# Uses 'corner' to make corner plots (cite dx.doi.org/10.5281/zenodo.45906)

import filefinder as ff
import numpy as np
import matplotlib.pyplot as plt

import itertools

dumpfile = raw_input("What is the dump filename? ")
nspiralmax = input("Maximum number of spirals to plot: ")    
# Get list of spiral files

filenames = ff.find_sorted_local_input_fileset(dumpfile+'_spiral*.dat')

# number of spirals = number of files

nspiral = len(filenames)

# Figure markers

markers = itertools.cycle((',', '+', '.', 'o','^'))
colours = itertools.cycle(('red','green', 'orange', 'blue', 'purple'))
# Set up figure

fig1=plt.figure()
ax1 = fig1.add_subplot(111)
ax1.set_xlabel('x')
ax1.set_ylabel('y')
#ax1.set_xlim(-40,40)
#ax1.set_ylim(-40,40)

ispiral = 0


# For each file in the list:
for filename in filenames:
        
    ispiral = ispiral + 1
    if (ispiral > nspiralmax): break
    print 'Reading input filename ',filename
    
    # Read in the spiral data
    data = np.genfromtxt(filename)
            
    xi = data[:,0]
    yi = data[:,1]
    zi = data[:,2]
        
    
    ax1.scatter(xi,yi, color=colours.next(), marker=markers.next())
    ax1.annotate(str(ispiral), xy=(xi[-1],yi[-1]), arrowprops=dict(facecolor='black', shrink=0.05), horizontalalignment='right', verticalalignment='bottom')
        

plt.show()    
    
