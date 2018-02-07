# Written 12/06/17 by dh4gan
# Reads in the output from FORTRAN 90 code 'spiralfind' after fitting with `fit_spiraldata.py'
# (List of spiral files, with x y z points)

# Plots the spiral data plus best MCMC fits 


import filefinder as ff
import numpy as np
import matplotlib.pyplot as plt
import io_spiral

import itertools

# Choose which spiral to fit

spiralchoice,spiraltext,nparams = io_spiral.choose_spiral()

if(spiralchoice =='logarithmic'):
    generate_curve = io_spiral.generate_logspiral_curvem  

elif(spiralchoice == 'varlogarithmic'):
    generate_curve = io_spiral.generate_varlogspiral_curvem

elif(spiralchoice=='hyperbolic'):
    generate_curve = io_spiral.generate_hypspiral_curvem    

elif(spiralchoice=='power'):
    generate_curve = io_spiral.generate_powspiral_curvem

elif(spiralchoice=='rpitch'):
    optfunc = io_spiral.generate_rpitchspiral_curvem

numcol = 1
xsigncol = nparams+3
ysigncol = xsigncol+1
nplot = 1000


# Choose files to plot
#dumpfile = raw_input("What is the dump filename? ")
nspiralmax = input("Maximum number of spirals to plot: ")    



# Load dumpfile names from spirallist.txt
dumpfiles = np.loadtxt('spirallist.txt', dtype='string',skiprows=1)

try:
    nfiles = len(dumpfiles)
    print "There are ",nfiles, " dumpfiles"
    print dumpfiles
except TypeError:
    print "There is only one dump to analyse: ",dumpfiles
    dumpfiles = np.array([dumpfiles, ' '])

for dumpfile in dumpfiles:

    # Get list of spiral files
    if dumpfile ==' ': continue

    filenames = ff.find_sorted_local_input_fileset(dumpfile+'_spiral*.dat')

    # number of spirals = number of files

    nspiral = len(filenames)

    # Figure markers
    
    markers = itertools.cycle((',', '+', '.', 'o','^'))
    colours = itertools.cycle(('red','green', 'orange', 'blue', 'purple'))

    # Set up figure

    fig1=plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.set_xlabel('x (AU)',fontsize=22)
    ax1.set_ylabel('y (AU)',fontsize=22)

    ispiral = 0

    # Read list of fit parameters for all spirals

    fitfile = dumpfile+'_spirals.chiminfits'

    print 'Reading fitted parameters from file ', fitfile
    
    fitdata = np.genfromtxt(fitfile,skiprows=2)

    print fitdata

    # For each file in the list:
    for filename in filenames:
        
        ispiral = ispiral + 1
        if (ispiral > nspiralmax): break
        print 'Reading input filename ',filename
    
        # Read in the spiral data
        data = np.genfromtxt(filename)
            
        if fitdata[ispiral-1,numcol]<2: 
            print 'Skipping: low spiral count'
            continue
    
            #    if fitdata[ispiral-1,chicol]>10.0:
            #        print 'Skipping: High chi^2'
            #        continue
            
        xi = data[:,0]
        yi = data[:,1]
        zi = data[:,2]

        mfit = fitdata[ispiral-1,2:numcol+1+nparams]
                            
        xsignfit = fitdata[ispiral-1,xsigncol]
        ysignfit = fitdata[ispiral-1,ysigncol]

        print mfit,xsignfit,ysignfit
        xplot,yplot = generate_curve(xi[1],yi[1],xi[-1],yi[-1],mfit,xsign=xsignfit,ysign=ysignfit,nplot=nplot)
    
        
        ax1.plot(xplot,yplot, color='red')
        ax1.scatter(xi,yi, color='green', marker='x')
        
    
    fig1.savefig(dumpfile+'_spirals_fitted.png')    
    
        

    
    
