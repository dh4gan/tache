# Given input parameters, plot example logarithmic spiral

import numpy as np
import matplotlib.pyplot as plt
from io_spiral import rpitchspiral_theta

npoints = 1000
a = 1.0
d1 = 0.1
d2 = 1.5
rp = 10
eta = 0.25

x0 = 0.0
y0 = 0.0

rmin = 1.0
rmax = 10.0*rp

r = np.linspace(rmin,rmax,num=npoints)

print "r dependent pitch spiral parameters:"
print "a: ",a
print "d1,d2 : ",d1,d2
print "rp: ", rp

fig1 = plt.figure()
ax1  = fig1.add_subplot(311)
ax2 = fig1.add_subplot(312)
ax3 = fig1.add_subplot(313)

x = np.zeros(npoints)
y = np.zeros(npoints)
theta = np.zeros(npoints)
pitch = np.zeros(npoints)
b = np.zeros(npoints)

for i in range(npoints):
	x[i],y[i],theta[i],pitch[i],b[i] = rpitchspiral_theta(r[i],a,d1,d2,eta,rp,x0,y0)

ax1.plot(x,y)
ax1.set_xlabel('x')
ax1.set_ylabel('y')


ax2.plot(r,theta)
ax2.set_xlabel('r')
ax2.set_ylabel('theta')

ax3.plot(r/rp,pitch*180.0/np.pi)
ax3.set_xlabel('r/rp')
ax3.set_ylabel('pitch')
ax3.set_xscale('log')
ax3.set_yscale('log')

plt.show()
