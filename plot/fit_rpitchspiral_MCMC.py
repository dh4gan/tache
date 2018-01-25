# Written 02/02/17 by dh4gan
# Reads in the output from'spiralfind'
# (List of spiral files, with x y z points)

# Runs an MCMC (Metropolis-Hastings) to find best fitting spiral with 
# r-dependent pitch angle (cf Zhu et al 2015)

import filefinder as ff
import numpy as np
import matplotlib.pyplot as plt
import sys
from io_spiral import opt_chisquared_rpitchspiral

import corner as c

npoints = 100 # number of evaluations to find point distance from spiral model
nburn = 1000 # Burn in sequence length
nanalyse = 10 # Only fit this many arms
nsamples = int(1.0e5) # Total number of MCMC samplings
nsubsample = 10 # Rate of subsampling of MCMC 

# R-dependent pitch spiral has 7  model parameters
#m = [a,hp,alpha,eta,rp,x0,y0]
nparams = 7

# Parameters initially sampled uniformly during burn in period (flat prior)

mmin = [1.0,   0.0, -5.0, 0.01, 300.0,-10.0,-10.0]
mmax = [100.0, 1.5,  5.0,  0.5, 500.0, 10.0, 10.0]

# After burn in, MC sampling from multivariate Gaussian

sigma_m = [0.01, 0.01, 0.01, 0.01, 1.0, 0.01,0.01]

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


    if dumpfile == ' ': continue
    # Get list of spiral files

    filenames = ff.find_sorted_local_input_fileset(dumpfile+'_spiral*.dat')

    ispiral = 0

    # For each file in the list:
    for filename in filenames:
        
        ispiral = ispiral + 1
	if ispiral > nanalyse: 
	    break
        print 'Reading input filename ',filename
    
        # Read in the spiral data
        data = np.genfromtxt(filename)
            
        xi = data[:,0]
        yi = data[:,1]
        zi = data[:,2]
        
        # Begin MCMC process

        mvalues = []
        chivalues = []
        xsignvalues = []
        ysignvalues = []
    
        isample = 0
 
        minit = [1.0, 1.0, 0.0, 400.0, 0.25, 0.0, 0.0]
        xsigninit= 1.0
        ysigninit = 1.0
    
        chimin = opt_chisquared_rpitchspiral(minit,xi, yi, npoints,xsign=xsigninit, ysign=ysigninit,sigma=0.1,verbose=False)    
    
        m = minit
        mnext = minit
        xsign = xsigninit
        ysign = ysigninit
        
        globalchimin = 1.0e30
        naccept = 0.0
        
        print 'Beginning MCMC spiral curve fitting '
        print 'Burning in: burn sequence length is ',nburn
        while isample < nsamples+nburn:
            
            # choose parameters for next sample
        
            # Explore the parameter space widely initially
            if(isample < nburn):  
                for k in range(nparams):
                    mnext[k] = (mmax[k]-mmin[k])*np.random.rand() + mmin[k]
            
            # After burn in, choose next variable from multivariate Gaussian
            else:
                for k in range(nparams):
                    mnext[k] = np.random.randn()*sigma_m[k] + m[k]
                    
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
    
            chinext = opt_chisquared_rpitchspiral(mnext,xi, yi, npoints, xsign=xsignnext, ysign=ysignnext,sigma=0.1,verbose=False)
                        
            # Calculate likelihood ratio
        
            ratio = np.exp(-chinext + chimin)
        
            # If this exceeds the test value, then accept point
            randtest = np.random.rand()        
            isample = isample + 1
        
                        
            if(randtest<ratio):
                            
                chimin = chinext
                m = mnext
#                a = anext
#                d1 = d1next
#		d2 = d2next
#                eta = etanext
#		rp = rpnext
#                x0 = x0next
#                y0= y0next
                xsign = xsignnext
                ysign = ysignnext
                
                if(chimin<globalchimin):
                    globalchimin = chimin
                    global_m_min = m
#                globalamin = a 
#                globald1min = d1
#		globald2min = d2
#                globaletamin = eta
#		globalrpmin = rp 
#                globalx0min = x0
#                globaly0min = y0
                globalxsign = xsign
                globalysign = ysign
            
                if(isample > nburn):
 
                    naccept = naccept + 1  

                    print 'Sample: %5i %+4.2e %+4.2e %+4.2e %+4.2e %+4.2e %+4.2e +%4.2e %i %i %4.2e' % (isample,m[0],m[1],m[2],m[3],m[4],m[5],m[6],m[7],xsign,ysign,chimin)  
                    mvalues.append(m)                
                    chivalues.append(chimin)
            
        # Save values to MCMC output file
    
        print naccept, ' samples were accepted out of ',(nsamples-nburn)
        rejection_rate = 100.0*(nsamples-nburn-naccept)/(nsamples-nburn)
        print 'Rejection rate is ',rejection_rate
        
        nsubsample = 10
        print 'Subsampling by a factor of ',nsubsample
    
        nMCMC = len(mvalues)                
        MCMCfile = filename+'.MCMC'
    
        mvalues = np.array(mvalues)
    
        # Subsample the total
 
        msubsample = mvalues[::nsubsample,:]
    
        np.savetxt(MCMCfile, msubsample)
        
        print 'MCMC minimum: '
        print global_m_min[:],globalxsign, globalysign
    
        # Write best fit to array for later file write

        spiralfit = global_m_min
        spiralfit.insert(0,len(xi))
        spiralfit.insert(0,ispiral)
        spiralfit.append(chimin)
        spiralfit.append(globalxsign)
        spiralfit.append(globalysign)

        MCMCfits.append(spiralfit)
    

    # end of loop over spirals

    # Write MCMC fit parameters to file
    outputfile = dumpfile+'_spirals.MCMCfits'
    outputformat = '%3i %3i'
    outputformat = outputformat +' %+7.5e '*(nparams+1)+'%i %i'
    
    np.savetxt(outputfile,MCMCfits,header="Best fitting parameters for log spiral via MCMC \n",fmt=outputformat)
    
    
