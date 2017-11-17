import numpy as np
import matplotlib.pyplot as plt
import read_eigenvalue_file as eigen
import read_eigenvector_file as eigvec


def read_eigenvalue_file(filename):
    '''Calls the f2py functions needed to read the eigenvalue file'''

    nelement = eigen.find_number_entries(filename)
    x,y,z,eigenelement,eigenvalues = eigen.read_file(filename,nelement)

    eigenvalues = eigenvalues.transpose()

    return nelement,x,y,z,eigenelement,eigenvalues

def read_eigenvector_file(filename):
    '''Calls the f2py functions needed to read the eigenvector file'''

    nelement = eigvec.find_number_entries(filename)
    eigenvecelement,eigenvectors = eigvec.read_file(filename,nelement)

    return nelement,eigenvecelement,eigenvectors

def find_corresponding_eigenvector_file(eigenvaluefile):
   '''Given an input eigenvalue file, finds corresponding eigenvector file'''
   
   # eigenvalue file format = "eigenvalues"+(T/V)+"_"+filename

   filename = eigenvaluefile[11:]

   eigenvectorfile = "eigenvectors"+filename
   return eigenvectorfile


def find_corresponding_eigenvalue_file(eigenvectorfile):
    '''Given an input eigenvector file, finds corresponding eigenvalue file'''
    
    filename = eigenvectorfile[12:]
    eigenvaluefile = "eigenvalues"+filename
    return eigenvaluefile


def classify_eigenvalue(eigenvalues, threshold):
    '''Given 3 eigenvalues, and some threshold, returns an integer
    'iclass' corresponding to the number of eigenvalues below the threshold
    
    iclass = 0 --> clusters (3 +ve eigenvalues, 0 -ve)
    iclass = 1 --> filaments (2 +ve eigenvalues, 1 -ve)
    iclass = 2 --> sheet (1 +ve eigenvalues, 2 -ve)
    iclass = 3 --> voids (0 +ve eigenvalues, 3 -ve)
    
    '''
    
    iclass = 0
    
    for i in range(3):

        if(eigenvalues[i]<threshold):
            iclass +=1
            
    return int(iclass)

def classify_all_eigenvalues(eigenvalues,nelement,threshold):
    '''Classifies `nelement` eigenvalues using classify_eigenvalue'''

    classification = np.empty(nelement,dtype='int')

    for i in range(nelement):
        eig = eigenvalues[i,:]
        classification[i] = classify_eigenvalue(eig,threshold)


    # Also return array indices for each class

    return classification

def print_class_counts(classification,filename):

     # Array indices for each class
    nclusters = np.size(classification[classification[:]==0])
    nfilaments = np.size(classification[classification[:]==1])
    nsheets = np.size(classification[classification[:]==2])
    nvoids = np.size(classification[classification[:]==3])

    print "File: ",filename
    print "Clusters: ",nclusters
    print "Filaments: ",nfilaments
    print "Sheets: ",nsheets
    print "Voids: ",nvoids

    return nclusters,nfilaments,nsheets,nvoids


def plot_all_classes_xy(x,y,filename,classification,threshold):

    xmin = np.amin(x)
    xmax = np.amax(x)

    ymin = np.amin(y)
    ymax = np.amax(y)

    # Array indices for each class
    index_clusters = classification[:]==0
    index_filaments = classification[:]==1
    index_sheets = classification[:]==2
    index_voids = classification[:]==3

    fig1, ((axcl,axfil),(axsh,axvd)) = plt.subplots(2,2, sharex='col',sharey='row')

    # Decide on plot title using filename:
    # Files containing "T" have been calculated using the tidal tensor
    # Files containing "V" have been calculated using the velocity shear tensor
    
    plot_title = "Tensor Classification, Threshold ="+str(threshold)
        
    if "T" in filename:
        plot_title = "Tidal Tensor Classification, Threshold ="+str(threshold)

    if "V" in filename:
        plot_title = "Velocity Shear Tensor Classification, Threshold ="+str(threshold)

    plt.suptitle(plot_title)
    axcl.set_title("Clusters")
    axcl.hexbin(x[index_clusters],y[index_clusters],gridsize=500)
    axcl.set_xlim(xmin,xmax)
    axcl.set_ylim(ymin,ymax)

    axfil.set_title("Filaments")
    axfil.hexbin(x[index_filaments],y[index_filaments],gridsize=500)
    axfil.set_xlim(xmin,xmax)
    axfil.set_ylim(ymin,ymax)
    
    axsh.set_title("Sheets")
    axsh.hexbin(x[index_sheets],y[index_sheets],gridsize=500)
    axsh.set_xlim(xmin,xmax)
    axsh.set_ylim(ymin,ymax)
    
    axvd.set_title("Voids")
    axvd.hexbin(x[index_voids],y[index_voids],gridsize=500)
    axvd.set_xlim(xmin,xmax)
    axvd.set_ylim(ymin,ymax)
    

    return fig1

    #outputfile = "classify_"+filename+"_"+str(threshold)+".png"
    #print "Plotting to file ",outputfile
         
    #plt.savefig(outputfile,format="png")
