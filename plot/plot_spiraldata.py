# Written 02/02/17 by dh4gan
# Reads in the output from FORTRAN 90 code 'spiralfind'
# (List of spiral files, with x y z points)

# Runs an MCMC (Metropolis-Hastings) to find best fitting logarithmic spiral
# Uses 'corner' to make corner plots (cite dx.doi.org/10.5281/zenodo.45906)

import filefinder as ff
import numpy as np
import matplotlib.pyplot as plt

from io_spiral import get_chisquared,logspiral_x,logspiral_y,find_minimum_t

import corner as c


npoints = 100 # number of evaluations to find point distance from spiral model
nburn = 50 # Burn in sequence length

nsamples = int(1.0e5) # Total number of MCMC samplings

# Parameters initially sampled uniformly during burn in period (flat prior)

amin = 3.5
amax = 5.5

bmin = 0.1
bmax = 0.3

x0min = -2.0
x0max = 2.0

y0min = -2.0
y0max = 2.0

# After burn in, MC sampling from multivariate Gaussian

sigma_a = 0.01
sigma_b = 0.01
sigma_x0 = 0.01
sigma_y0 = 0.01
    
    
#dumpfiles = ['clus300', 'clus301', 'clus302','clus303']
dumpfiles = ['clus300']

for dumpfile in dumpfiles:
    # Get list of spiral files

    filenames = ff.find_sorted_local_input_fileset(dumpfile+'_spiral*.dat')

    # Set up file containing all fit parameters

    f_fit = open(dumpfile+'_spirals.fits','w')
    line = '# Fits for individual spirals'

    # Set up final plot (data + fits)

    fig1=plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.set_xlabel('x (AU)')
    ax1.set_ylabel('y (AU)')

    ispiral = 0

    # For each file in the list:
    for filename in filenames:
        
        ispiral = ispiral + 1
        print 'Reading input filename ',filename
    
        # Read in the spiral data
        data = np.genfromtxt(filename)
            
        xi = data[:,0]
        yi = data[:,1]
        zi = data[:,2]
        
        # Begin MCMC process
    
        avalues = []
        bvalues = []
        x0values = []
        y0values = []
        chivalues = []
        xsignvalues = []
        ysignvalues = []
    
        isample = 0
    
        ainit = 1.0
        binit = 1.0
        x0init = 0.0
        y0init = 0.0
        xsigninit= 1.0
        ysigninit = 1.0
    
        chimin = get_chisquared(xi, yi, ainit, binit, x0init, y0init, npoints,xsign=xsigninit, ysign=ysigninit,sigma=0.1)    
    
        a=ainit
        b=binit
        x0=x0init
        y0=y0init
        xsign = xsigninit
        ysign = ysigninit
        
        globalchimin = 1.0e30
        naccept = 0.0
        
        print 'Beginning MCMC spiral curve fitting '
        print 'Burning in: burn sequence length is ',nburn
        while isample < nsamples:
            
            # choose parameters for next sample
        
            # Explore the parameter space widely initially
            if(isample < nburn):                
                anext = (amax-amin)*np.random.rand() + amin
                bnext = (bmax-bmin)*np.random.rand() + bmin
                x0next = (x0max-x0min)*np.random.rand() + x0min
                y0next = (y0max-y0min)*np.random.rand() + y0min
            
            # After burn in, choose next variable from multivariate Gaussian
            else:
                anext = np.random.randn()*sigma_a + a
                bnext = np.random.randn()*sigma_b + b
                x0next = np.random.randn()*sigma_x0 + x0
                y0next = np.random.randn()*sigma_y0 + y0
                    
            xsignnext = 1.0
            randtest = np.random.rand()
            if(randtest>0.5):
                xsignnext = -1.0
             
            ysignnext = 1.0
            randtest = np.random.rand()
            if(randtest>0.5):
                ysignnext = -1.0
    
            #xsignnext = 1.0
            #ysignnext = -1.0
    
            chinext = get_chisquared(xi, yi, anext, bnext, x0next, y0next, npoints, xsign=xsignnext, ysign=ysignnext,sigma=0.5)
                        
            # Calculate likelihood ratio
        
            ratio = np.exp(-chinext + chimin)
        
            # If this exceeds the test value, then accept point
            randtest = np.random.rand()        
            isample = isample + 1
        
                        
            if(randtest<ratio):
                            
                chimin = chinext
                a = anext
                b = bnext
                x0 = x0next
                y0= y0next
                xsign = xsignnext
                ysign = ysignnext
                
            if(chimin<globalchimin):
                globalchimin = chimin
                globalamin = a 
                globalbmin = b 
                globalx0min = x0
                globaly0min = y0
                globalxsign = xsign
                globalysign = ysign
            
            if(isample > nburn):
 
                naccept = naccept + 1           
                print 'Sample: ', isample,a,b,x0,y0,xsign,ysign,chimin
                avalues.append(a)
                bvalues.append(b)
                x0values.append(x0)
                y0values.append(y0)
                xsignvalues.append(xsign)
                ysignvalues.append(ysign)
                chivalues.append(chimin)
            
        # Save values to MCMC output file
    
        print naccept, ' samples were accepted out of ',(nsamples-nburn)
        rejection_rate = 100.0*(nsamples-nburn-naccept)/(nsamples-nburn)
        print 'Rejection rate is ',rejection_rate
        
        nsubsample = 10
        print 'Subsampling by a factor of ',nsubsample
    
        nMCMC = len(avalues)                
        MCMCfile = filename+'.MCMC'
    
        allsamples = np.array((avalues,bvalues,x0values,y0values))    
        allsamples = np.transpose(allsamples)
    
        # Subsample the total
    
        print allsamples.shape
        allsamples = allsamples[::nsubsample,:]
    
    
        print allsamples.shape
    
        np.savetxt(MCMCfile, allsamples)
    
        #f_obj = open(MCMCfile,'w')
    
        #for i in range(nMCMC):
        #    line = str(avalues[i])+ '  '+str(bvalues[i]) + '   '+str(x0values[i]) + '   '+str(y0values[i]) + '   '+str(chivalues[i])
        #    line = line + '   '+str(xsignvalues[i])+'    '+str(ysignvalues[i])+ ' \n'
        #    f_obj.write(line)
        
        #f_obj.close()
    
    
        print 'MCMC minimum: '
        print globalamin, globalbmin, globalx0min, globaly0min, globalxsign,globalysign
    
        # Write best fit to file
    
        line = str(ispiral)+'   '+str(len(xi))+ '  '+ str(globalamin)+ '  '+str(globalbmin) + '   '+str(globalx0min) + '   '+str(globaly0min) + '   '+str(chimin)
        line = line + '   '+str(globalxsign)+'    '+str(globalysign)+ ' \n'
    
        f_fit.write(line)
    
        # Create corner plot for MCMC parameters
    
        cplot = c.corner(allsamples, labels=['$a$','$b$','$x_0$ (AU)','$y_0$ (AU)'], label_kwargs = {"fontsize": 22})
    
        cplot.savefig(dumpfile+'_MCMC_corner_'+str(ispiral)+'.png')
        
        # Find minimum and maximum t for spiral plotting
    
        tmin,sepmin = find_minimum_t(xi[0],yi[0], globalamin,globalbmin,globalx0min,globaly0min, npoints,xsign=globalxsign,ysign=globalysign)
        tmax,sepmin = find_minimum_t(xi[-1],yi[-1], globalamin,globalbmin,globalx0min,globaly0min, npoints,xsign=globalxsign,ysign=globalysign)
        
        print tmin, tmax
        nplot = 100
        t = np.linspace(tmin,tmax,num=nplot)
        
        xplot = np.zeros(nplot)
        yplot = np.zeros(nplot)
    
        for i in range(nplot):
        
            xplot[i] = logspiral_x(t[i], globalamin, globalbmin, globalx0min,xsign=globalxsign)
            yplot[i] = logspiral_y(t[i], globalamin, globalbmin, globaly0min,ysign=globalysign)
        
        
        ax1.plot(xplot,yplot, color='red')
        ax1.scatter(xi,yi, color='green', marker='x')
    
    f_fit.close()
    fig1.savefig(dumpfile+'_spirals_fitted.png')    
    
