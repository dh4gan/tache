# Loads MCMC output from fit_rpitchspiral and makes a corner plot

import numpy as np
import corner as c
import filefinder as ff

nfiles = input("How many files to generate corner plots for?")

MCMCfiles = ff.find_sorted_local_input_fileset('*.dat.MCMC')

print "Generating ",nfiles, " corner plots (for logarithmic fits)"

for i in range(nfiles):
    MCMCfile = MCMCfiles[i]
    allsamples = np.genfromtxt(MCMCfile)

    # Create corner plot for MCMC parameters

    cplot = c.corner(allsamples, labels=['$a$','$\delta_1$','$\alpha$','$\eta$','$r_p$','$x_0$ (AU)','$y_0$ (AU)'], label_kwargs = {"fontsize": 22})

    cplot.savefig(MCMCfile+'.rpitch.cornerplot.png')
