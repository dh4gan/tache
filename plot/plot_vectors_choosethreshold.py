# Written 9/10/14 by dh4gan
# Code reads in  output "eigen" file from sph_tidal_tensor
# Applies a user-defined threshold to the particles' tensor eigenvalues
# to classify the region they inhabit
# Classification is plotted in colour-coded format

import numpy as np
import matplotlib.pyplot as plt
import io_tache as io
import filefinder as ff
from scipy.interpolate import griddata
from scipy.stats import maxwell,rayleigh


def safeArcCos(x):
    '''A safe version of arccos that handles abs(x)>1'''
    
    x[x>1.0] = 1.0
    x[x<-1.0] = -1.0    
        
    return np.arccos(x)

# Column indices for different properties
xcol = 0
ycol = 1
zcol = 2
eigcols = range(3,6)

# Colours of different classes
# Black= Void
# Green = filament
# Blue = sheet
# Red = cluster

colourchoice = [[0,0,0],[0,1,0],[0,0,1],[1,0,0]]


vectorfile = ff.find_local_input_files('eigenvectors*')
valuefile = io.find_corresponding_eigenvalue_file(vectorfile)
threshold = input("What is the threshold for classification? ")

# Read in eigenvalue file

print "Reading eigenvalue file ", valuefile

nelement,x,y,z,eigenpart,eigenvalues = io.read_eigenvalue_file(valuefile)


# Read in eigenvector file

print "Reading eigenvector file ", vectorfile

npart, eigenvecpart,eigenvectors = io.read_eigenvector_file(vectorfile)

if(not np.array_equal(eigenpart,eigenvecpart)): 
    print "WARNING: files have inconsistent particle IDs"

tryagain = 'y'

# Now enter a do while loop until user is happy with result

while(tryagain!='n'):

    print "Classifying eigenvalues according to threshold ", threshold  
 
    classification = io.classify_all_eigenvalues(eigenvalues,npart,threshold)

    # Give some statistics on the classification

    nclusters,nfilaments,nsheets,nvoids = io.print_class_counts(classification,valuefile)

    xfil = []
    yfil=[]
    zfil = []
    filamentvectors = []
    
    xsheet = []
    ysheet = []
    zsheet = []
    sheetnormals = []
    
    print "Determining filament flows and sheet normals"

    for i in range(npart):
        
        if classification[i]==1:
            # If particle is classified as a filament, then 
            # negative eigenvalue's associated
            # eigenvector should give flow direction
        
            ifilament = np.argmin(eigenvalues[i,:])
        
            eigvec = eigenvectors[ifilament,:,i]                    
            
            
            filamentvectors.append(eigvec)
            xfil.append(x[i])
            yfil.append(y[i])
            zfil.append(z[i])

        if classification[i]==2:
            # If particle is classified as a sheet, 
            # then +ve eigenvalue's associated
            # eigenvector should give sheet normal
        
            isheet = np.argmax(eigenvalues[i,:])
        
            eigvec = eigenvectors[isheet,:,i]
                        
            sheetnormals.append(eigvec)
            xsheet.append(x[i])
            ysheet.append(y[i])
            zsheet.append(z[i])     
         
    # Now plot this data

    xmin = np.amin(x)
    xmax = np.amax(x)

    ymin = np.amin(y)
    ymax = np.amax(y)
    
    zmin = np.amin(z)
    zmax = np.amax(z)

    #xmin = -300.0
    #xmax = 300.0

    #ymin = -300.0
    #ymax = 300.0
    
    # Interpolate the vector fields
    
    xx = np.arange(xmin,xmax,0.1)
    yy = np.arange(xmin,xmax,0.1)
    zz = np.arange(zmin,zmax,0.1)
        
    xgrid,ygrid,zgrid =np.meshgrid(xx,yy,zz)

    points = np.column_stack((xfil,yfil,zfil))
    filamentvectors = np.array(filamentvectors)
    
    xvalues = filamentvectors[:,0]
    yvalues = filamentvectors[:,1]
    zvalues = filamentvectors[:,2]
            
    print 'Interpolating Filament eigenvector Field in 3D'        
    
    gridxfil = griddata(points,xvalues,(xgrid,ygrid,zgrid),method='linear',fill_value=0.0)
    gridyfil = griddata(points,yvalues,(xgrid,ygrid,zgrid),method='linear',fill_value=0.0)
    gridzfil = griddata(points,zvalues,(xgrid,ygrid,zgrid),method='linear',fill_value=0.0)
    
    nx = gridxfil.shape[0]
    ny = gridxfil.shape[1]
    nz = gridxfil.shape[2]            
    
    print 'Normalising Vectors'
    
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                magnitude = np.sqrt(gridxfil[i,j,k]*gridxfil[i,j,k] + gridyfil[i,j,k]*gridyfil[i,j,k] + gridzfil[i,j,k]*gridzfil[i,j,k])
                if(magnitude > 0.0):
                    gridxfil[i,j,k] = gridxfil[i,j,k]/magnitude
                    gridyfil[i,j,k] = gridyfil[i,j,k]/magnitude
                    gridzfil[i,j,k] = gridzfil[i,j,k]/magnitude
                else:
                    gridxfil[i,j,k] = np.nan
                    gridyfil[i,j,k] = np.nan
                    gridzfil[i,j,k] = np.nan
    
    
    points = np.column_stack((xsheet,ysheet,zsheet))
    sheetnormals = np.array(sheetnormals)
    
    xvalues = sheetnormals[:,0]
    yvalues = sheetnormals[:,1]
    zvalues = sheetnormals[:,2]
    
    print 'Interpolating Sheet eigenvector Field in 3D'
    
    gridxsheet = griddata(points,xvalues,(xgrid,ygrid,zgrid),method='linear',fill_value=0.0)
    gridysheet = griddata(points,yvalues,(xgrid,ygrid,zgrid),method='linear',fill_value=0.0)
    gridzsheet = griddata(points,zvalues,(xgrid,ygrid,zgrid),method='linear',fill_value=0.0)
    
    print 'Normalising Vectors'
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                magnitude = np.sqrt(gridxsheet[i,j,k]*gridxsheet[i,j,k] + gridysheet[i,j,k]*gridysheet[i,j,k] + gridzsheet[i,j,k]*gridzsheet[i,j,k])
                if(magnitude > 0.0):
                    gridxsheet[i,j,k] = gridxsheet[i,j,k]/magnitude
                    gridysheet[i,j,k] = gridysheet[i,j,k]/magnitude
                    gridzsheet[i,j,k] = gridzsheet[i,j,k]/magnitude
                else:
                    gridxsheet[i,j,k] = np.nan
                    gridysheet[i,j,k] = np.nan
                    gridzsheet[i,j,k] = np.nan
            
    print 'Done'
            
    # Calculate angles between the vectors (dot product)
    
    print 'Calculating angles between sheet and filament vectors'
        
    angles = np.zeros(gridxfil.shape,dtype=float)    
    
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):            
                                                
                angles[i,j,k] = gridxfil[i,j,k]*gridxsheet[i,j,k]+gridyfil[i,j,k]*gridysheet[i,j,k] + gridzfil[i,j,k]*gridzsheet[i,j,k]
                #print i,j, gridxfil[i,j], gridxsheet[i,j], gridyfil[i,j],gridysheet[i,j], angles[i,j]
        
    #angles = angles[np.logical_not(np.isnan(angles))]
    angles = safeArcCos(angles)
    
    # Attempt to fit angle data to Maxwell distribution
    
    print 'Fitting to Maxwell and Rayleigh distributions'
    
    maxwell_params = maxwell.fit(angles[np.logical_not(np.isnan(angles))])
    rayleigh_params = rayleigh.fit(angles[np.logical_not(np.isnan(angles))])
    #rice_params = rice.fit(0.0,angles[np.logical_not(np.isnan(angles))])
    
    x = np.linspace(0.0,np.pi, 100)
    
    maxwell_curve = maxwell.pdf(x,loc=maxwell_params[0], scale=maxwell_params[1])
    rayleigh_curve = rayleigh.pdf(x,loc=rayleigh_params[0],scale=rayleigh_params[1])    
    
    #angles = np.nan_to_num(angles)
    
    
    print 'Done - Fit Parameters: ' 
    print 'Maxwell Distribution: ',maxwell_params    
    print 'Rayleigh Distribution: ',rayleigh_params    
            
    # Plot Vector Fields
    
    print 'Plotting Vector Fields'
    fig1, (axfil,axsh) = plt.subplots(2,sharex=True)

    markersize = 10.0
    
    # Decide on plot title using filename:
    # Files containing "T" have been calculated using the tidal tensor
    # Files containing "V" have been calculated using the velocity shear tensor
    
    plot_title = "Tensor Classification, Threshold ="+str(threshold)
    
    if "T" in valuefile:
        plot_title = "Tidal Tensor Classification, Threshold ="+str(threshold)
    
    if "V" in valuefile:
        plot_title = "Velocity Shear Tensor Classification, Threshold ="+str(threshold)

    plt.suptitle(plot_title)

    #plotxmin = -3.099
    #plotxmax = -1.875
    
    #plotymin = -5.486
    #plotymax = -4.268

    plotxmin =xmin
    plotxmax = xmax
    plotymin = ymin
    plotymax = ymax

    axfil.set_title("Filament Vectors")
    axfil.set_xlabel('x (pc)')
    axfil.set_ylabel('y (pc)')
    axfil.streamplot(xgrid[:,:,nz/2],ygrid[:,:,nz/2],gridxfil[:,:,nz/2],gridyfil[:,:,nz/2],density=1, arrowsize=2, color='blue')
    axfil.set_xlim(plotxmin,plotxmax)
    axfil.set_ylim(plotymin,plotymax)

    axsh.set_title("Sheet Normals")
    axsh.set_xlabel('x (pc)')
    axsh.set_ylabel('y (pc)')
    axsh.streamplot(xgrid[:,:,nz/2],ygrid[:,:,nz/2],gridxsheet[:,:,nz/2],gridysheet[:,:,nz/2], density=1, arrowsize=2,color='green')
    axsh.set_xlim(plotxmin,plotxmax)
    axsh.set_ylim(plotymin,plotymax)


    print 'Plotting Angles as a Scalar Field'

    fig2 = plt.figure()
    
    ax2 = fig2.add_subplot(111)    
    ax2.set_xlabel('x (pc)')
    ax2.set_ylabel('y (pc)')
    plt.pcolor(xgrid[:,:,nz/2],ygrid[:,:,nz/2],angles[:,:,nz/2], vmin = 0.0, vmax = np.pi)
    colourbar = plt.colorbar()
    colourbar.set_label('Angle between filament/sheet vectors (rad)')

    print 'Plotting distribution of angles'
    
    angles = angles.flatten()

    print np.amin(angles), np.amax(angles)
    fig3 = plt.figure()
    ax3 = fig3.add_subplot(111)
    ax3.set_ylabel('Relative Frequency')
    ax3.set_xlabel('Angle between Filament Vector and Sheet Normal (rad)')
    n,bins,patches = ax3.hist(angles[np.logical_not(np.isnan(angles))], bins=100, histtype='step', normed=True, label='Data')
    #ax3.plot(x,maxwell_curve, color='red',linewidth = 2, label='Maxwell Fit')
    #ax3.plot(x,rayleigh_curve, color='green',linewidth = 2, label='Rayleigh Fit')
    #ax3.legend()
    print np.amax(n), np.argmax(n), bins[np.argmax(n)]
    
    plt.show()
    
    # Check if user wants to try again
    
    tryagain = raw_input("Do you want to try a different threshold (y/n)?")
    
    if(tryagain!='n'):
        threshold = input("Enter a new threshold: ")
        
# End While Loop

# Once user is happy with the threshold, output this data to file

outputfile1 = "streamlines_"+vectorfile+"_"+str(threshold)+".png"
outputfile2 = "anglefield_"+vectorfile+"_"+str(threshold)+".png"
outputfile3 = "anglehist_"+vectorfile+"_"+str(threshold)+".png"

print 'Saving streamlines to file ', outputfile1
fig1.savefig(outputfile1,format="png")

print 'Saving angle field to file ', outputfile2
fig2.savefig(outputfile2,format="png")

print 'Saving angle histogram to file ', outputfile3
fig3.savefig(outputfile3,format="png")
