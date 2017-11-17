# Written 1/10/17 by dh4gan
# Code reads in output "eigenvalues" file from tache
# Applies a user-defined threshold to the particles' tensor eigenvalues
# to classify the region they inhabit, then plots them by azimuth
# Do this to attempt to fit spiral arm widths

import numpy as np
import matplotlib.pyplot as plt
import io_tache as io
from scipy.optimize import curve_fit


def func(x, *params):
    y = np.zeros_like(x)
    
    for i in range(0,len(params),3):
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]        
        
        y = y + amp*np.exp(-((x-ctr)**2/(2.0*(wid)**2)))
    return y


colourchoice = [[0,0,0],[0,1,0],[0,0,1],[1,0,0]]
classlabel = ['Cluster', 'Arm', 'Interarm', 'Void']

# Read threshold value from command line or as argument

filename = ff.find_local_input_files('eigenvalues*')
threshold = input("What is the threshold for classification? ")

print "Reading eigenvalue file ", filename

npart,x,y,z,eigenpart,eigenvalues = io.read_eigenvalue_file.(filename)


print "Classifying eigenvalues according to threshold ", threshold    
classification = io.classify_all_eigenvalues(eigenvalues,npart,threshold)

# Give some statistics on the classification

nclusters,nfilaments,nsheets,nvoids = io.print_class_counts(classification,filename)


# compute azimuths
rpart = np.zeros(npart)
phipart = np.zeros(npart)

for i in range(npart):
    eig = eigenvalues[i,:]
    rpart[i] = np.sqrt(x[i]*x[i]+y[i]*y[i])
    phipart[i] = np.arctan2(y[i],x[i])
    

    
# Now plot this data

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

rmin = 19.5
rmax = 20.0

# Remove particles outside of the radius of interest
classification[rpart[:]<rmin] = -10
classification[rpart[:]>rmax] = -10
    
ax1.set_xlabel(r'$\theta$ (rad) ', fontsize=20)
ax1.set_ylabel(r'Particle Count ', fontsize=20)    
for iclass in range(1,3):
    n, bins,patches = ax1.hist(phipart[classification[:]==iclass], histtype = 'step',bins=50, label=classlabel[iclass],linewidth=2)
        
    # Attempt to fit arm widths
    if(iclass==1):
            
        guess = [-2.0, 1000, 0.1, -0.75, 1000,0.1, 1,1000,0.1, 2.5,1000,0.1]
           
        bin_centers = bins[:-1] + 0.5*(bins[1:]-bins[:-1])
        popt, pcurve = curve_fit(func, bin_centers, n, p0=guess)
        
        print bin_centers, n
        print popt
            
        fit = func(bin_centers, *popt)
                        
        ax1.plot(bin_centers, fit, label='fit to arms', color='red',linestyle='dashed')
        
        
        #ax1.hist(rpart[classification[:]==iclass], alpha=0.1)
    
handles, labels = ax1.get_legend_handles_labels()
    
handles = [handles[2],handles[1],handles[0]]
labels = [labels[2],labels[1],labels[0]]
ax1.legend(handles,labels)    
    
plt.show()
fig1.savefig('class_vs_theta.png')
    
  
