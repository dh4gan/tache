# Written 12/06/17 by dh4gan
# Reads in the output from FORTRAN 90 code 'spiralfind' after fitting with `fit_spiraldata.py'
# (List of spiral files, with x y z points)

# Plots the spiral data plus best MCMC fits 


import filefinder as ff
import numpy as np
import matplotlib.pyplot as plt
from io_spiral import find_minimum_t_logspiral,logspiral_x,logspiral_y

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
        
    afit = fitdata[ispiral-1,acol]
    bfit = fitdata[ispiral-1,bcol]
    xofit = fitdata[ispiral-1,xocol]
    yofit = fitdata[ispiral-1,yocol]
    xsignfit = fitdata[ispiral-1,xsigncol]
    ysignfit = fitdata[ispiral-1,ysigncol]

    
    # Find minimum and maximum t for spiral plotting
    
    tmin,sepmin1 = find_minimum_t_logspiral(xi[2],yi[2], afit,bfit,xofit,yofit,npoints,xsign=xsignfit,ysign=ysignfit)
    tmax,sepmin2 = find_minimum_t_logspiral(xi[-1],yi[-1],afit,bfit,xofit,yofit,npoints,xsign=xsignfit,ysign=ysignfit)
        
                
    print tmin, tmax, sepmin1,sepmin2
    
#    if (tmax-tmin) > 5.0:
#        print 'Skipping: loopy fit'
#        continue
    
    nplot = 100
    t = np.linspace(tmin,tmax,num=nplot)
        
    xplot = np.zeros(nplot)
    yplot = np.zeros(nplot)
    
    for i in range(nplot):
        
        xplot[i] = logspiral_x(t[i], afit,bfit,xofit,xsign=xsignfit)
        yplot[i] = logspiral_y(t[i], afit,bfit,yofit,ysign=ysignfit)
        
        
    ax1.plot(xplot,yplot, color='red')
    ax1.scatter(xi,yi, color='green', marker='x')
    
#plt.show()    
fig1.savefig(dumpfile+'_spirals_fitted.png')    
    
        

    
    
