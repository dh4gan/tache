# Written 02/02/17 by dh4gan
# Reads in the output from FORTRAN 90 code 'spiralfind'
# (List of spiral files, with x y z points)

# Uses scipy.optimise to find the minimum chisquared to find best fitting spiral

import filefinder as ff
import numpy as np
import sys
import io_spiral
import scipy.optimize

#import corner as c

# Give initial guess parameters
# m = [a,b,x0,y0]

# Load dumpfile names from spirallist.txt
dumpfiles = np.loadtxt('spirallist.txt', dtype='string',skiprows=1)

nanalyse = 15

# Choose which spiral to fit

spiralchoice,spiraltext,nparams = io_spiral.choose_spiral()

# Also pick function to minimise from


if(spiralchoice =='logarithmic'):
    optfunc = io_spiral.opt_chisquared_logspiral
    minit = [10.0,0.1,0,0]

elif(spiralchoice=='hyperbolic'):
    optfunc = io_spiral.opt_chisquared_hypspiral
    minit = [10.0,0.0,0.0]

elif(spiralchoice=='power'):
    optfunc = io_spiral.opt_chisquared_powspiral
    minit = [10.0,1.0,0.0,0.0]

elif(spiralchoice=='rpitch'):
    optfunc = io_spiral.opt_chisquared_rpitchspiral
    minit = [150.0, 1.0, 2.0, 1.0, 400.0, 0.0, 0.0]



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

    # Set up array to hold fit parameters
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

        m = minit

        npoints = 100 # Number of evaluations to find point distance in spiral model 
        xsign = 1.0
        ysign = -1.0


        mopt = scipy.optimize.minimize(optfunc,m,args=(xi,yi,npoints,xsign,ysign),method='Nelder-Mead') 

        print mopt

        spiralfits = mopt.x.tolist()
        spiralfits.insert(0,len(xi))
        spiralfits.insert(0,ispiral)
        spiralfits.append(mopt.fun)
        spiralfits.append(xsign)
        spiralfits.append(ysign)
	chiminfits.append(spiralfits)

    print chiminfits
    outputfile = dumpfile+'_spirals.chiminfits'
    outputformat = '%3i %3i'
    outputformat = outputformat +' %+7.5e '*(nparams+1)+'%i %i'


    np.savetxt(outputfile,chiminfits,header="Minimum chisquared fits for "+spiraltext+" \n",fmt=outputformat)

print "Done"                    

    
    
