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

elif(spiralchoice=='hyperbolic'):
    generate_curve = io_spiral.generate_hypspiral_curvem    

elif(spiralchoice=='power'):
    generate_curve = io_spiral.generate_powspiral_curvem

elif(spiralchoice=='rpitch'):
    optfunc = io_spiral.generate_rpitchspiral_curvem

# Choose files to plot
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
ax1.set_xlabel('x (AU)',fontsize=22)
ax1.set_ylabel('y (AU)',fontsize=22)
#ax1.set_xlim(-120,120)
#ax1.set_ylim(-120,120)

#ax1.set_xlabel('x (kpc)',fontsize=22)
#ax1.set_ylabel('y (kpc)',fontsize=22)
#ax1.set_xlim(-5,5)
#ax1.set_ylim(-5,5)


ispiral = 0


# Read list of fits

numcol = 1
acol = 2
bcol = 3
xocol = 4
yocol = 5
chicol = 6
xsigncol = 7
ysigncol = 8

fitfile = dumpfile+'_spirals.chiminfits'
npoints = 100

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

    print mfit
    xplot,yplot = generate_curve(xi[2],yi[2],xi[-1],yi[-1],mfit,xsign=xsignfit,ysign=ysignfit,npoints=npoints)
    
        
    ax1.plot(xplot,yplot, color='red')
    ax1.scatter(xi,yi, color='green', marker='x')
    
#plt.show()    
fig1.savefig(dumpfile+'_spirals_fitted.png')    
    
        

    
    
