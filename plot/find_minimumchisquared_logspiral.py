# Written 02/02/17 by dh4gan
# Reads in the output from FORTRAN 90 code 'spiralfind'
# (List of spiral files, with x y z points)

# Uses scipy.optimise to find the minimum chisquared to find best fitting logarithmic spiral
# Uses 'corner' to make corner plots (cite dx.doi.org/10.5281/zenodo.45906)

import filefinder as ff
import numpy as np
#import matplotlib.pyplot as plt
import sys
from io_spiral import opt_chisquared_logspiral,logspiral_x,logspiral_y,find_minimum_t_logspiral
import scipy.optimize

#import corner as c

# Give initial guess parameters
# m = [a,b,x0,y0]

# Load dumpfile names from spirallist.txt
dumpfiles = np.loadtxt('spirallist.txt', dtype='string',skiprows=1)

nanalyse = 3

try:
    nfiles = len(dumpfiles)
    print "There are ",nfiles, " dumpfiles"
    print dumpfiles
except TypeError:
    print "There is only one dump to analyse: ",dumpfiles
    dumpfiles = np.array([dumpfiles, ' '])

for dumpfile in dumpfiles:

    if dumpfile == ' ': continue
    # Get list of spiral files

    filenames = ff.find_sorted_local_input_fileset(dumpfile+'_spiral*.dat')

    # Set up file containing all fit parameters

    f_fit = open(dumpfile+'_spirals.chimin','w')
    line = '# Minimised chisquared for individual spirals'

    chiminfits=[]
    ispiral = 0

    for filename in filenames:
	ispiral = ispiral+1
	if(ispiral>nanalyse):
	    break
	print "Minimising chi squared for file ", filename
        # Read in the spiral data
        data = np.genfromtxt(filename)

        xi = data[:,0]
        yi = data[:,1]

        m = np.zeros(4)
        m[0] = 20.0
        m[1] = 0.15
        m[2] = 0.0
        m[3] = 0.0    

        npoints = 100 # Number of evaluations to find point distance in spiral model 
        xsign = 1.0
        ysign = -1.0


        mopt = scipy.optimize.minimize(opt_chisquared_logspiral,m,args=(xi,yi,npoints,xsign,ysign),method='Nelder-Mead') 

        print mopt

	chiminfits.append([ispiral,len(xi),mopt.x[0],mopt.x[1],mopt.x[2],mopt.x[3],mopt.fun,xsign,ysign])


    print chiminfits
    outputfile = dumpfile+'_spirals.chimin'
    np.savetxt(outputfile,chiminfits,header="Minimum chisquared fits for log spirals \n")

print "Done"                    

    
    
