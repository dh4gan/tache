# Given input parameters, plot example logarithmic spiral

import numpy as np
import matplotlib.pyplot as plt
from io_spiral import logspiral_x, logspiral_y

npoints = 100
a = 84.0
pitch = 7.9


pitch = pitch* np.pi/180.0
b = 1.0/np.tan(0.5*np.pi-pitch)

x0 = 0.0
y0 = 0.0

tmin = 0.0
tmax = 10.0

t = np.linspace(tmin,tmax,num=npoints)


print "logarithmic spiral parameters:"
print "a: ",a
print "b (pitch angle): ",b, "(",pitch*180.0/np.pi,")"

fig1 = plt.figure()
ax1  = fig1.add_subplot(111)

x = np.zeros(npoints)
y = np.zeros(npoints)

for i in range(npoints):
	x[i] = logspiral_x(t[i],a,b,x0)
	y[i] = logspiral_y(t[i],a,b,y0)

plt.plot(x,y)

plt.show()
