# Created by D. Forgan (02/05/2013) 
# Uses glob to find files for user input
# (Potentially saves user some typing)

import glob
import numpy as np
import re

def sort_nicely(l):
    """
    Sort the given list in a 'Natural' order
    (Ned Batchelder's Compact Human Python Sort)
    """
    
    convert = lambda text:int(text) if text.isdigit() else text
    
    alphanum_key = lambda key:[ convert(c) for c in re.split('([0-9]+)',key)]
    
    l.sort(key=alphanum_key)

def find_local_input_files(stringmatch):
    '''Given a matching string (e.g. '*.txt') the function will create a list
    of matches for the user to select (or type in an alternative)'''
    
    print 'Searching local directory for input files'
    filechoices = glob.glob(stringmatch)
    
    # Number of matches
    nmatch = len(filechoices)
    
    print 'Detected ', nmatch, ' potential inputfiles in this directory'
    print 'Here are the options: '
    for i in range (nmatch):
        print '(',i+1,'): ', filechoices[i]
        
    if(nmatch>0): print 'If none of these files suit:'
    print '(',nmatch+1,'):  Manually enter a filename'
    
    userselect = input('Make a selection: ')
    
    if userselect==nmatch+1:
        filename = raw_input('Manually enter filename: ')
    else:
        filename = filechoices[userselect-1]
        
    print 'File ',filename, ' selected for read-in'
    return filename 


def find_local_input_fileset(stringmatch):
    '''Given a matching string (e.g. '*.txt') the function will create a list
    of all matches '''
    
    print 'Searching local directory for input files'
    filechoices = glob.glob(stringmatch)
    
    # Number of matches
    nmatch = len(filechoices)
    
    print 'Detected ', nmatch, ' potential inputfiles in this directory'    
    for i in range (nmatch):
        print '(',i+1,'): ', filechoices[i]
        
    return filechoices

def find_sorted_local_input_fileset(stringmatch):
    '''Given a matching string (e.g. '*.txt') the function will create a sorted list
    of all matches '''
    
    filechoices = glob.glob(stringmatch)
    
    # Number of matches
    nmatch = len(filechoices)
    
    # Sort using Natural sort (see function at top)
    sort_nicely(filechoices)
    for i in range (nmatch):
        print '(',i+1,'): ', filechoices[i]
    
    
    print 'Detected ', nmatch, ' potential inputfiles in this directory'    
    
    
    return filechoices
    

def decide_local_output_file(inputfile, stringmatch,fileformat='ps'):
    '''Given a matching string, the function gives a list of existing
    output files to write over, or allows an alternative to be entered)'''
    
    print 'Searching local directory for possible output filenames'
    filechoices = glob.glob(stringmatch)
    
    print 'Creating default outputfile name'
    defaultfile = inputfile.rsplit('.',1)[0]
    defaultfile = defaultfile +'.'+fileformat
    
    filechoices.insert(0,defaultfile)
    
    # Number of matches
    nmatch = len(filechoices)
    
    print 'Detected ', nmatch, ' potential output files to overwrite in this directory'
    print 'Here are the options: '
    for i in range (nmatch):
        if(i==0): 
            print 'Default (',i+1,'): ', filechoices[i]
        else: 
            print '(',i+1,'): ', filechoices[i]
        
    if(nmatch>0): print 'If you wish to create a new file: '
    print '(',nmatch+1,'):  Manually enter a filename'
    
    userselect = input('Make a selection: ')
    
    if userselect==nmatch+1:
        filename = raw_input('Manually enter filename: ')
    else:
        filename = filechoices[userselect-1]
        
    print 'File ',filename, ' selected for output'
    return filename 




    
def find_local_numbered_input_files(ignores):
    '''Searches for lists of files numbered .1 etc
    ignores array forces the code to trim away extra text'''
    
    print 'Searching local directory for sets of numbered input filenames'
    
    numtry = 0
    filechoices = []
    print ignores, '*.'+str(numtry)+'*'
    # Search from '*.0*' to '*.9*' 
    
    while filechoices == []:
        numtry = numtry + 1
        filechoices = glob.glob('*.'+str(numtry)+'*')
    
    # Find set of unique file prefixes - strip numbered suffixes
    
    prefixes = []
    for i in  range(len(filechoices)):
        filechoices[i] = filechoices[i].rsplit('.',1)[0]
        
        
    filechoices = np.unique(filechoices)
    
    print filechoices
    
    for filename in filechoices:
        print 'Before stripping ',filename
        for phrase in ignores:                     
            if filename.endswith(phrase):        
                filename = filename[:len(filename)-len(phrase)]               
            
            print filename, phrase
        prefixes.append(filename)
    
    prefixes = np.unique(prefixes)
    nmatch = len(prefixes)
    
    
    # Now find maximum and minimum filenumbers for each prefix
    
    maxnum = []
    minnum = []
    filenums = []
    
    for pre in prefixes:
        filechoices = glob.glob(pre+'*')
    
        for i in range(len(filechoices)):
            num = filechoices[i].rsplit('.',1)[1]
            if num.isdigit():
                filenums.append(int(num))
            
        filenums = np.unique(filenums)        
                
        maxnum.append(np.amax(filenums))
        minnum.append(np.amin(filenums))
        filenums = []
        
    print 'Detected ', nmatch, ' set of potential input files in this directory'
    print 'Here are the options (prefix, min filenumber,max filenumber): '
    for i in range (nmatch):
        print '(',i+1,') : ', prefixes[i], minnum[i],maxnum[i]
        
    if(nmatch>0): print 'If you wish to create a new file: '
    print '(',nmatch+1,'):  Manually enter a filename'
    
    userselect = input('Make a selection: ')
    
    if userselect==nmatch+1:
        prefixchoose = raw_input('Manually enter filename: ')
    else:
        prefixchoose = prefixes[userselect-1]
        
    print 'File prefix',prefixchoose, ' selected for output'
    return prefixchoose 
