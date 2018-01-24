# Written 15/1/18 by dh4gan
# Script reads spiralmembers.dat file from spiralfind
# Also reads best fit parameters for each arm
# Then computes distance of particle from arm as a function of radius

import filefinder as ff
import io_tache
import io_spiral
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

npoints = 5000

print 'Select membership file to analyse:'

memberfile = ff.find_local_input_files('*_spiralmembers.dat')
ispiral = input("Which spiral to plot? ")

# Determine eigenvalue file name from memberfile name
eigenfile = memberfile[:-18]

# Read spiralmembers file
print 'Reading spiral membership in file',memberfile

x,y,z,spiralmember = io_spiral.read_spiralmembership(memberfile)
print 'Read ', len(spiralmember), ' elements'
#print spiralmember

# Read best fits (either .chimin or .fitparams)
fitparamfile = eigenfile+'_spirals.chimin'

fitdata = np.genfromtxt(fitparamfile,skiprows=2)

# Find fit parameters for selected arm

a = fitdata[ispiral-1,2]
b = fitdata[ispiral-1,3]
x0 = fitdata[ispiral-1,4]
y0 = fitdata[ispiral-1,5]
xsign = fitdata[ispiral-1,7]
ysign = fitdata[ispiral-1,8]

# Find all elements belonging to that arm

imember = spiralmember[:]==ispiral

x = x[imember]
y = y[imember]
z = z[imember]

xorigin = 0.0
yorigin = 0.0

nmember = len(x)

print 'Found ', nmember, ' members of spiral ', ispiral

# For each element:
# compute r, sepmin (minimum distance from spiral)
# save to arrays

#nmember = 1000
r = np.zeros(nmember)
t = np.zeros(nmember)
sepmin = np.zeros(nmember)
weight = np.zeros(nmember)

for i in range(nmember):

    r[i] = io_spiral.separation(xorigin,yorigin,x[i],y[i])
    t[i], sepmin[i] = io_spiral.find_minimum_t_logspiral(x[i],y[i],a,b,x0,y0,npoints,xsign=xsign,ysign=ysign)

    print i,r[i],t[i], sepmin[i]

tmin = np.amin(t)
tmax = np.amax(t)

weight[:] = 1.0/float(nmember)

print 'Minimum, maximum r: ', np.amin(r), np.amax(r)
print 'Minimum, maximum t: ', tmin, tmax
print 'Generating curve: '
print a,b,x0,y0,xsign,ysign


xspiral, yspiral = io_spiral.generate_logspiral_curve(tmin,tmax,a,b,x0,y0,xsign=xsign,ysign=ysign,npoints=npoints)


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
#ax2 = fig1.add_subplot(212)
ax1.set_xlabel('R (kpc)',fontsize=22)
ax1.set_ylabel('Spine Distance (kpc)',fontsize=22)
counts, xbins,ybins, image = ax1.hist2d(r/10.0,sepmin/10.0,bins=20, range=[[1.0,4.0],[0.0,0.1]],normed=False,cmap='rainbow')
#plt.colorbar(image,ax=ax1)

maxcount = counts.max()

print maxcount, np.median(counts)

clevels = [50,70,90,95]
clabels = [str(i)+'%' for i in clevels]
clevels = [np.percentile(counts[np.nonzero(counts)],i) for i in clevels]


print clevels
print clabels

CS= ax1.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],colors='white',levels=clevels)

fmt={}
for l,s in zip(CS.levels,clabels):
    fmt[l]=s

plt.clabel(CS,fontsize=16,fmt=fmt)

#ax1.hist(sepmin[:100])
#ax2.scatter(x[:100],y[:100])
#ax2.plot(xspiral,yspiral,color='red')


fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.set_ylabel('Relative Frequency',fontsize=22)
ax2.set_xlabel('Spine Distance (kpc)',fontsize=22)
ax2.hist(sepmin/10.0,bins=50, histtype='step', label='all radii', linewidth=2,normed=True)

sepclose  = sepmin[np.logical_and(r[:]>=10.0,r[:]<20.0)]
sepfar = sepmin[np.logical_and(r[:]>=20.0,r[:]<30.0)]
ax2.hist(sepclose/10.0, histtype = 'step',label = '$1.0 < r < 2.0 $ kpc',linewidth=2,normed=True)
ax2.hist(sepfar/10.0,histtype = 'step',label = '$2.0 < r < 3.0 $ kpc',linewidth=2,normed=True)

ax2.legend(loc='upper right')
plt.show()
fig1.savefig(eigenfile+'spiral_'+str(ispiral)+'width_vs_r.png')
fig2.savefig(eigenfile+'spiral_'+str(ispiral)+'width1D.png')

