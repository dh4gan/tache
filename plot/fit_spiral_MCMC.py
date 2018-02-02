# Written 02/02/17 by dh4gan
# Reads in output from 'spiralfind'
# (List of spiral files, with x y z points)

# Runs an MCMC (Metropolis-Hastings) to find best fitting logarithmic spiral (in x,y)
# Outputs samples and a file containing best fit parameters for each spiral

import filefinder as ff
import numpy as np
import sys
import io_spiral

npoints = 100 # number of evaluations to find point distance from spiral model
nburn = 500 # Burn in sequence length
nanalyse = 3 # Only analyse this many arms
nsamples = int(1.0e4) # Total number of MCMC samplings
nsubsample = 10 # Subsample rate of MCMC sampling


# Choose which spiral to fit
spiralchoice,spiraltext,nparams = io_spiral.choose_spiral()

# Must define parameter value ranges 
# Also define sigma for multivariate Gaussian for sampling

if(spiralchoice =='logarithmic'):
    opt_func = io_spiral.opt_chisquared_logspiral
    minit = [10.0,0.1,0.0,0.0]

    mmin = [1.0,   0.1, -10.0, -10.0]
    mmax = [100.0, 0.7,  10.0, 10.0]

    sigma_m = [0.01,0.01,0.01,0.01]

elif(spiralchoice=='hyperbolic'):
    opt_func = io_spiral.opt_chisquared_hypspiral
    minit = [10.0,0,0]

    mmin = [1.0,   -100.0, -100.0]
    mmax = [1.0e4,  100.0, 100.0]

    sigma_m = [0.1,0.1,0.1]


elif(spiralchoice=='rpitch'):
    opt_func = io_spiral.opt_chisquared_rpitchspiral
    minit = [150.0, 1.0, 2.0, 1.0, 400.0, 0.0, 0.0]

    mmin = [1.0,   0.0, -5.0, 0.01, 300.0,-10.0,-10.0]
    mmax = [100.0, 1.5,  5.0,  0.5, 500.0, 10.0, 10.0]

    sigma_m = [0.01, 0.01, 0.01, 0.01, 1.0, 0.01,0.01]


# Load dumpfile names from spirallist.txt
try:
    dumpfiles = np.loadtxt('spirallist.txt', dtype='string',skiprows=1)
except IOError:
    print "File spirallist.txt not found: quitting"
    exit()

try:
    nfiles = len(dumpfiles)
    print "There are ",nfiles, " dumpfiles"
    print dumpfiles
except TypeError:
    print "There is only one dump to analyse: ",dumpfiles
    dumpfiles = np.array([dumpfiles, ' '])

for dumpfile in dumpfiles:

    if dumpfile == ' ': continue
    # Get list of spiral files for this dump

    filenames = ff.find_sorted_local_input_fileset(dumpfile+'_spiral*.dat')

    MCMCfits = []
    ispiral = 0

    # For each spiral file in the list:
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
            
        xsigninit= 1.0
        ysigninit = 1.0
    
        chimin = opt_func(minit,xi, yi, npoints,xsign=xsigninit, ysign=ysigninit,sigma=0.1,verbose=False)    
    
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
        
            chinext = opt_func(mnext,xi, yi, npoints, xsign=xsignnext, ysign=ysignnext,sigma=0.1,verbose=False)
                        
            # Calculate likelihood ratio
        
            ratio = np.exp(-chinext + chimin)
        
            # If this exceeds the test value, then accept point
            randtest = np.random.rand()        
            isample = isample + 1
            print 'Sample: ','%5i'%isample,' '.join(['%4.2e']*len(m))%tuple(m), '%4.2e %i %i'%(chimin,xsign,ysign)
            if(randtest<ratio):
                            
                chimin = chinext
                m = mnext
                
                if(chimin<globalchimin):
                    globalchimin = chimin
                    global_m_min = m
                    globalxsign = xsign
                    globalysign = ysign
            
                if(isample > nburn):
                    naccept = naccept + 1         
                    mvalues.append(m)
            
        # Save values to MCMC output file
    
        print naccept, ' samples were accepted out of ',(nsamples-nburn)

        rejection_rate = 100.0*(nsamples-nburn-naccept)/(nsamples-nburn)
        print 'Rejection rate is ',rejection_rate
        
        print 'Subsampling by a factor of ',nsubsample
    
        nMCMC = len(mvalues)                
        MCMCfile = filename+'.MCMC'
    
        # Subsample the total
        mvalues = np.array(mvalues)

        msubsample = mvalues[::nsubsample,:]
            
        np.savetxt(MCMCfile, msubsample)
            
        print 'MCMC minimum: '
        print global_m_min[:], globalxsign,globalysign
    
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
    
    np.savetxt(outputfile,MCMCfits,header="Best fitting parameters for "+spiralxtext +" via MCMC \n",fmt=outputformat)
    
    
