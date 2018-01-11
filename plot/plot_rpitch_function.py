# Plot pitch angle function
import numpy as np
import matplotlib.pyplot as plt

from io_spiral import rpitch_phi

npoints = 100
rp = 450.0
a = 20
d1 = 10*np.pi/180.0
d2 = 1.5

print d1,d2, np.power(rp,d2)

r = np.linspace(0,2.0*rp, num=npoints) 
pitch = np.zeros(npoints)
b = np.zeros(npoints)

for i in range(npoints):
    pitch[i], b[i] = rpitch_phi(a,d1,d2,r[i],rp)



fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(r/rp,pitch)
ax1.set_xlabel(r'$(r/r_p)$', fontsize=22)
ax1.set_ylabel(r'$\phi(r)$', fontsize=22)

plt.show()
