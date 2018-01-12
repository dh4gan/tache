# Written 12/06/17 by dh4gan
# Reads in the output from FORTRAN 90 code 'spiralfind' after fitting with `fit_spiraldata.py'
# (List of spiral files, with x y z points)

# Plots the spiral data plus best MCMC fits 


import filefinder as ff
import numpy as np
import matplotlib.pyplot as plt
from io_spiral import find_minimum_t_rpitchspiral,rpitchspiral_theta

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
d1col = 3
d2col = 4
rpcol = 5
xocol = 4
yocol = 5
chicol = 6
xsigncol = 7
ysigncol = 8

fitfile = dumpfile+'_spirals.fits'
npoints = 100

print 'Reading fitted parameters from file ', fitfile
    
fitdata = np.genfromtxt(fitfile)

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
        
    afit = fitdata[ispiral-1,acol]
    d1fit = fitdata[ispiral-1,d1col]
    d2fit = fitdata[ispiral-1,d2col]
    rpfit = fitdata[ispiral-1,rpcol]
    xofit = fitdata[ispiral-1,xocol]
    yofit = fitdata[ispiral-1,yocol]
    xsignfit = fitdata[ispiral-1,xsigncol]
    ysignfit = fitdata[ispiral-1,ysigncol]

    
    # Find minimum and maximum r for spiral plotting
    
    rmin = np.sqrt(xi[2]*xi[2] + yi[2]*yi[2])
    rmax = np.sqrt(xi[-1]*xi[-1] + yi[-1]*yi[-1])
    
#    if (tmax-tmin) > 5.0:
#        print 'Skipping: loopy fit'
#        continue
    
    nplot = 100
    rplot = np.linspace(rmin,rmax,num=nplot)
        
    xplot = np.zeros(nplot)
    yplot = np.zeros(nplot)
    thetaplot[i] = np.zeros(nplot)
    for i in range(nplot):
        
        xplot[i],yplot[i],thetaplot[i] = rpitchspiral_x(rplot[i], afit,d1fit,d2fit,rpfit,xofit,yofit,xsign=xsignfit,ysign=ysignfit)        
        
    ax1.plot(xplot,yplot, color='red')
    ax1.scatter(xi,yi, color='green', marker='x')
    
#plt.show()    
fig1.savefig(dumpfile+'_rpitchspirals_fitted.png')    
    
        

    
    
