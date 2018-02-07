import numpy as np

#
# Functions for I/O etc
#


spiralchoices = ['logarithmic','varlogarithmic','hyperbolic','power','rpitch']
spiraltexts = ['Logarithmic Spiral', 'Logarithmic Spiral, varying pitch','Hyperbolic Spiral','Power Spiral', 'r-dependent pitch spiral']
nspiralchoices = len(spiralchoices)
nspiralparams = [4,5,3,4,7]

def choose_spiral():

    userselect = nspiralchoices+5
    print 'Choose which type of spiral to analyse: '

    while(userselect > nspiralchoices):

        print 'Here are the options: ',nspiralchoices
        for i in range (nspiralchoices):
            print '(',i+1,'): ', spiralchoices[i]
            
        userselect = input('Make a selection: ')
    
    if userselect>nspiralchoices:
        print "Choice out of range: please try again"
    

    spiralchoice = spiralchoices[userselect-1]
    spiraltext = spiraltexts[userselect-1]
    nparams = nspiralparams[userselect-1]
    print spiraltext, " has been selected"
    return spiralchoice,spiraltext,nparams


def separation(x1,y1,x2,y2):
    '''Return separation of x y coordinates'''
    return np.sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))

def read_spiralmembership(filename):
    '''Reads the spiral membership data of all particles analysed'''
    spiraldata = np.genfromtxt(filename)
    x = spiraldata[:,1]
    y = spiraldata[:,2]
    z = spiraldata[:,3]
    spiralmember = spiraldata[:,4]
    return x,y,z,spiralmember


############################################
# Functions for fitting logarithmic spirals
############################################


# Functions to deliver parametric form of spiral (x)
def logspiral_x(t,a,b,x0,xsign=1):    
    return xsign*a*np.exp(b*t)*np.cos(t) + x0

# Overloaded function so that a single array of model parameters can be passed
def logspiral_xm(t,m,xsign=1):
    m[0] = a
    m[1] = b
    m[2] = x0

    return logspiral_x(t,a,b,x0,xsign=xsign)

# Same for y
def logspiral_y(t,a,b,y0,ysign=1):    
    return ysign*a*np.exp(b*t)*np.sin(t) + y0

def logspiral_ym(t,m,ysign=1):
    m[0] = a
    m[1] = b
    m[3] = y0

    return logspiral_y(t,a,b,y0,ysign=ysign)
    

# Find the minimum t value for (xi,yi) given spiral parameters (a,b,x0,y0)
# Multiple local minima possible, must be careful

def find_minimum_t_logspiral(xi,yi, a,b,x0,y0, npoints,xsign=1.0,ysign=1.0):
    '''Find the minimum t value for points (xi,yi) given spiral parameters (a,b,x0,y0)'''
    
    t = np.linspace(0.0,10.0,num=npoints)
    
    tmin = -1.0
    sepmin = 1.0e30
    
    for i in range(npoints):
        
        x = logspiral_x(t[i], a, b, x0, xsign=xsign)
        y = logspiral_y(t[i], a, b, y0, ysign=ysign)
        
        sep = separation(xi, yi, x, y)                
        #print t[i],sep,x,y,sepmin
        if(sep<sepmin):
            sepmin = sep
            tmin = t[i]
        

    return tmin,sepmin


#
# Generate x,y points for a logarithmic spiral curve
#
def generate_logspiral_curve(xbegin,ybegin,xend,yend,a,b,x0,y0,xsign=1,ysign=1,nplot=100):
    '''Generate x,y, points of a logarithmic spiral'''

    nfind = 100
    tmin, sepmin1 = find_minimum_t_logspiral(xbegin,ybegin,a,b,x0,y0,nfind,xsign=xsign,ysign=ysign)

    tmax, sepmin2 = find_minimum_t_logspiral(xend,yend,a,b,x0,y0,nfind,xsign=xsign,ysign=ysign)

    #print tmin, xbegin,ybegin, sepmin1
    #print tmax, xend, yend, sepmin2

   # print a,b,x0,y0
    t = np.linspace(tmin,tmax,num=nplot)
    xspiral = np.zeros(nplot)
    yspiral = np.zeros(nplot)

    for i in range(nplot):
        xspiral[i] = logspiral_x(t[i],a,b,x0,xsign=xsign)
        yspiral[i] = logspiral_y(t[i],a,b,y0,ysign=ysign)
        #print t[i], xspiral[i],yspiral[i]
    
    return xspiral,yspiral

def generate_logspiral_curvem(xbegin,xend,ybegin,yend,m,xsign=1,ysign=1,nplot=100):
    '''Generate x,y, points of a logarithmic spiral'''
    a = m[0]
    b = m[1]
    x0 = m[2]
    y0 = m[3]
    return generate_logspiral_curve(xbegin,xend,ybegin,yend,a,b,x0,y0,xsign,ysign,nplot)




def get_chisquared_logspiral(x,y,a,b,x0,y0,npoints,xsign=1.0,ysign=1.0,sigma=1.0):
    '''Returns the chi-squared of a logarithmic spiral model given arrays x,y
    Assumes uniform errors'''
    
    tmin = np.zeros(len(x))
    sepmin = np.zeros(len(x))
        
    for i in range(len(x)):     
        tmin[i], sepmin[i] = find_minimum_t_logspiral(x[i], y[i], a, b, x0, y0, npoints,xsign=xsign,ysign=ysign)    
    
    return np.sum(sepmin)/(2.0*len(x)*sigma*sigma)


def opt_chisquared_logspiral(m,x,y,npoints,xsign=1.0,ysign=1.0,sigma=1.0,verbose=True):
    '''Wrapper for scipy.optimize: the chi-squared of a logarithmic spiral model given arrays x,y
    Assumes uniform errors'''
    
    a = m[0]
    b = m[1]
    x0 = m[2]
    y0 = m[3]

    chisquared = get_chisquared_logspiral(x,y,a,b,x0,y0,npoints,xsign,ysign,sigma)

    if(verbose):print 'Chi^2 - %4.2f'%chimin,': Parameters  ',' '.join(['%4.2e']*len(m))%tuple(m)
    return chisquared

#
# End of functions for logarithmic spirals
#

##############################################
# Functions for log spiral with varying pitch
##############################################
def varpitch_phi(a,b0,n,r):

    #b = b0 *pow(r/a,n)
    b = b0 + n*(r-a)
    pitch = np.arctan(b)
    return pitch,b

def varlogspiral_theta(r,a,b0,n,x0,y0,xsign=1,ysign=1):
    pitch, b = varpitch_phi(a,b0,n,r) 
#    b =0.3    
    r0 = np.sqrt(x0*x0 + y0*y0)
    theta = np.log((r-r0)/a)/b
    theta = np.mod(theta,2.0*np.pi)

    x = xsign*r*np.cos(theta)
    y = ysign*r*np.sin(theta)
    return x,y,theta,pitch,b

def varlogspiral_thetam(r,m,xsign=1,ysign=1):

    a = m[0]
    b0 = m[1]
    n = m[2] 
    x0 = m[3]
    y0 = m[4]
    
    
    return varlogspiral_theta(r,a,b0,n,x0,y0,xsign=xsign,ysign=ysign)

def varlogspiral_x(t,a,b0,n,r,x0,xsign=1):
	# Compute b parameter
        pitch,b = varpitch_phi(a,b0,n,r)

	# Now compute spiral x position
	return xsign*a*np.exp(b*t)*np.cos(t) + x0

def varlogspiral_xm(t,m,r,xsign=1):
    a = m[0] 
    b0 = m[1]
    n = m[2] 
    x0 = m[3]
    y0 = m[4]    
    
    return varlogspiral_x(t,a,b0,n,r,x0,xsign)

def varlogspiral_y(t,a,b0,n,r,y0,ysign=1):
	# Compute b parameter
        pitch,b = varpitch_phi(a,b0,n,r)

	# Now compute spiral x position
	return ysign*a*np.exp(b*t)*np.sin(t) + y0

def varlogspiral_ym(t,m,r,ysign=1):
    a = m[0] 
    b0 = m[1]
    n = m[2] 
    y0 = m[4]    
    
    return varlogspiral_y(t,a,b0,n,r,y0,ysign)


def find_minimum_t_varlogspiral(xi,yi, a,b0,n,x0,y0, npoints,xsign=1.0,ysign=1.0):
    '''Find the minimum t value for points (xi,yi) given spiral parameters (a,b,x0,y0)'''

    t = np.linspace(0.0,10.0,num=npoints)

    tmin = -1.0
    sepmin = 1.0e30

    for i in range(npoints):

      
	r = np.sqrt(xi*xi + yi*yi)
        x = varlogspiral_x(t[i], a, b0,n,r,x0, xsign=xsign)
        y = varlogspiral_y(t[i], a, b0,n,r,y0,ysign=ysign)

        sep = separation(xi, yi, x, y)

        if(sep<sepmin):
            sepmin = sep
            tmin = t[i]

    return tmin,sepmin


#
# Generate x,y points for an rpitch  spiral curve
#
def generate_varlogspiral_curve(xbegin,ybegin,xend,yend,a,b0,n,x0,y0,xsign=1,ysign=1,nplot=100):
    '''Generate x,y, points of an rpitch spiral'''

    nfind = 1000

    rmin = np.sqrt(xbegin*xbegin + ybegin*ybegin)
    rmax = np.sqrt(xend*xend + yend*yend)
    print rmin, rmax
    r = np.linspace(rmin,rmax,num=nplot)
    xspiral = np.zeros(nplot)
    yspiral = np.zeros(nplot)
    theta = np.zeros(nplot)
    pitch = np.zeros(nplot)
    b = np.zeros(nplot)

    for i in range(nplot):

        xspiral[i], yspiral[i],theta,pitch,b = varlogspiral_theta(r[i],a,b0,n,x0,y0,xsign=xsign,ysign=ysign)
    
    return xspiral,yspiral

def generate_varlogspiral_curvem(xbegin,xend,ybegin,yend,m,xsign=1,ysign=1,nplot=100):
    '''Generate x,y, points of an rpitch spiral'''
    a = m[0] 
    b0 = m[1]
    n = m[2] 
    x0 = m[3]
    y0 = m[4]    
  
    return generate_varlogspiral_curve(xbegin,xend,ybegin,yend,a,b0,n,x0,y0,xsign,ysign,nplot)



def get_chisquared_varlogspiral(x,y,a,b0,n,x0,y0,npoints,xsign=1.0,ysign=1.0,sigma=1.0):

    '''Returns the chi-squared of a r-dependent pitch spiral model given arrays x,y
    Assumes uniform errors'''

    tmin = np.zeros(len(x))
    sepmin = np.zeros(len(x))

    for i in range(len(x)):        
        tmin[i], sepmin[i] = find_minimum_t_varlogspiral(x[i], y[i], a, b0,n,x0, y0, npoints,xsign=xsign,ysign=ysign)

    return np.sum(sepmin)/(2.0*len(x)*sigma*sigma)

    
def opt_chisquared_varlogspiral(m,x,y,npoints,xsign=1.0,ysign=1.0,sigma=1.0,verbose=True):

    a = m[0]
    b0 = m[1]
    n = m[2]
    x0 = m[3]
    y0 = m[4]

    chimin = get_chisquared_varlogspiral(x,y,a,b0,n,x0,y0,npoints,xsign,ysign,sigma)

    if(verbose):print 'Chi^2 - %4.2f'%chimin,': Parameters  ',' '.join(['%4.2e']*len(m))%tuple(m)
    return chimin



##########################################
# Functions for hyperbolic spiral
#########################################

# Parametric Functions

def hypspiral_x(t,c,x0,xsign=1):    
    return xsign*c*np.cos(t)/t + x0

def hypspiral_y(t,c,y0,ysign=1):    
    return ysign*c*np.sin(t)/t + y0

def hypspiral_xm(t,m,xsign=1):    
    c = m[0]
    x0 = m[1]
    return hyperbolic_spiralx(t,c,x0,xsign)

def hypspiral_ym(t,m,ysign=1):    
    c = m[0]
    y0 = m[2]
    return hyperbolic_spiralx(t,c,y0,ysign)


#
# Generate x,y points for a hyperbolic spiral curve
#
def generate_hypspiral_curve(xbegin,ybegin,xend,yend,c,x0,y0,xsign=1,ysign=1,nplot=100):
    '''Generate x,y, points of a hyperbolic spiral'''

    nfind = 1000
    tmin, sepmin1 = find_minimum_t_hypspiral(xbegin,ybegin,c,x0,y0,nfind,xsign=xsign,ysign=ysign)

    nfind = 1000
    tmax, sepmin2 = find_minimum_t_hypspiral(xend,yend,c,x0,y0,nfind,xsign=xsign,ysign=ysign)

    print tmin, tmax
    print sepmin1,sepmin2

    t = np.linspace(tmin,tmax,num=nplot)
    xspiral = np.zeros(nplot)
    yspiral = np.zeros(nplot)

    for i in range(nplot):
        xspiral[i] = hypspiral_x(t[i],c,x0,xsign=xsign)
        yspiral[i] = hypspiral_y(t[i],c,y0,ysign=ysign)
    
    return xspiral,yspiral

def generate_hypspiral_curvem(xbegin,xend,ybegin,yend,m,xsign=1,ysign=1,nplot=100):
    '''Generate x,y, points of a hyperbolic spiral'''
    c = m[0]
    x0 = m[1]
    y0 = m[2]

    return generate_hypspiral_curve(xbegin,xend,ybegin,yend,c,x0,y0,xsign,ysign,nplot)


# Find the minimum t value for (xi,yi) given spiral parameters (a,b,x0,y0)
# Multiple local minima possible, must be careful

def find_minimum_t_hypspiral(xi,yi, c,x0,y0, npoints,xsign=1.0,ysign=1.0):
    '''Find the minimum t value for points (xi,yi) given spiral parameters (a,b,x0,y0)'''
    
    t = np.linspace(0.0,10000.0,num=npoints)
    
    tmin = -1.0
    sepmin = 1.0e30
    
    for i in range(npoints):
        
        x = hypspiral_x(t[i], c, x0, xsign=xsign)
        y = hypspiral_y(t[i], c, y0, ysign=ysign)
        
        sep = separation(xi, yi, x, y)                
        
        if(sep<sepmin):
            sepmin = sep
            tmin = t[i]
        

    return tmin,sepmin

def get_chisquared_hypspiral(x,y,c,x0,y0,npoints,xsign=1.0,ysign=1.0,sigma=1.0):
    '''Returns the chi-squared of a hyperbolic spiral model given arrays x,y
    Assumes uniform errors'''
    
    tmin = np.zeros(len(x))
    sepmin = np.zeros(len(x))
        
    for i in range(len(x)):     
        tmin[i], sepmin[i] = find_minimum_t_hypspiral(x[i], y[i], c, x0, y0, npoints,xsign=xsign,ysign=ysign)    
    
    return np.sum(sepmin)/(2.0*len(x)*sigma*sigma)


def opt_chisquared_hypspiral(m,x,y,npoints,xsign=1.0,ysign=1.0,sigma=1.0,verbose=True):
    '''Wrapper for scipy.optimize: the chi-squared of a hyperbolic spiral model given arrays x,y
    Assumes uniform errors'''
    
    c = m[0]    
    x0 = m[1]
    y0 = m[2]

    chisquared = get_chisquared_hypspiral(x,y,c,x0,y0,npoints,xsign,ysign,sigma)
    

    if(verbose):print 'Chi^2 - %4.2f'%chimin,': Parameters  ',' '.join(['%4.2e']*len(m))%tuple(m)

    return chisquared

#########################################
# End of hyperbolic spiral functions
#########################################


############################################################################
# Power spiral functions (r = a*theta^n)
# (n=1: Archimedes spiral)
# (n=1/2: Fermat Spiral)
############################################################################

# Functions to deliver parametric form of spiral (x)
def powspiral_x(t,a,n,x0,xsign=1):    
    return xsign*a*pow(t,n)*np.cos(t) + x0

# Overloaded function so that a single array of model parameters can be passed
def powspiral_xm(t,m,xsign=1):
    a = m[0]
    n = m[1]
    x0 = m[2]

    return powspiral_x(t,a,n,x0,xsign=xsign)

# Same for y
def powspiral_y(t,a,n,y0,ysign=1):    
    return ysign*a*pow(t,n)*np.sin(t) + y0

def powspiral_ym(t,m,ysign=1):
    a = m[0]
    n = m[1]
    y0 = m[3]

    return powspiral_y(t,a,n,y0,ysign=ysign)
    
#
# Generate x,y points for a power spiral curve
#
def generate_powspiral_curve(xbegin,ybegin,xend,yend,a,n,x0,y0,xsign=1,ysign=1,nplot=100):
    '''Generate x,y, points of a power spiral'''

    nfind = 1000
    tmin, sepmin1 = find_minimum_t_powspiral(xbegin,ybegin,a,n,x0,y0,nfind,xsign=xsign,ysign=ysign)
    tmax, sepmin2 = find_minimum_t_powspiral(xend,yend,a,n,x0,y0,nfind,xsign=xsign,ysign=ysign)

    print tmin, tmax
    t = np.linspace(tmin,tmax,num=nplot)
    xspiral = np.zeros(nplot)
    yspiral = np.zeros(nplot)

    for i in range(nplot):
        xspiral[i] = powspiral_x(t[i],a,n,x0,xsign=xsign)
        yspiral[i] = powspiral_y(t[i],a,n,y0,ysign=ysign)
    
    return xspiral,yspiral

def generate_powspiral_curvem(xbegin,xend,ybegin,yend,m,xsign=1,ysign=1,nplot=100):
    '''Generate x,y, points of a logarithmic spiral'''
    a = m[0]
    n = m[1]
    x0 = m[2]
    y0 = m[3]
    return generate_powspiral_curve(xbegin,xend,ybegin,yend,a,n,x0,y0,xsign,ysign,nplot)


# Find the minimum t value for (xi,yi) given spiral parameters (a,b,x0,y0)
# Multiple local minima possible, must be careful

def find_minimum_t_powspiral(xi,yi, a,n,x0,y0, npoints,xsign=1.0,ysign=1.0):
    '''Find the minimum t value for points (xi,yi) given spiral parameters (a,b,x0,y0)'''
    
    t = np.linspace(-1.0,10.0,num=npoints)
    
    tmin = -1.0
    sepmin = 1.0e30
    
    for i in range(npoints):
        
        x = powspiral_x(t[i], a, n, x0, xsign=xsign)
        y = powspiral_y(t[i], a, n, y0, ysign=ysign)
        
        sep = separation(xi, yi, x, y)                
        
        if(sep<sepmin):
            sepmin = sep
            tmin = t[i]
        

    return tmin,sepmin

def get_chisquared_powspiral(x,y,a,n,x0,y0,npoints,xsign=1.0,ysign=1.0,sigma=1.0):
    '''Returns the chi-squared of a power spiral model given arrays x,y
    Assumes uniform errors'''
    
    tmin = np.zeros(len(x))
    sepmin = np.zeros(len(x))
        
    for i in range(len(x)):     
        tmin[i], sepmin[i] = find_minimum_t_powspiral(x[i], y[i], a, n, x0, y0, npoints,xsign=xsign,ysign=ysign)    
    
    return np.sum(sepmin)/(2.0*len(x)*sigma*sigma)


def opt_chisquared_powspiral(m,x,y,npoints,xsign=1.0,ysign=1.0,sigma=1.0,verbose=True):
    '''Wrapper for scipy.optimize: the chi-squared of a logarithmic spiral model given arrays x,y
    Assumes uniform errors'''
    
    a = m[0]
    n = m[1]
    x0 = m[2]
    y0 = m[3]

    chimin = get_chisquared_powspiral(x,y,a,n,x0,y0,npoints,xsign,ysign,sigma)

    if(verbose):print 'Chi^2 - %4.2f'%chimin,': Parameters  ',' '.join(['%4.2e']*len(m))%tuple(m)

    return chimin

#
# End of functions for power spirals
#


############################################################################
# Functions for spiral with r-dependent pitch angle ('rpitch')
# Expected for tidally driven arms by low mass companions in low mass discs
# (Zhu et al (2015), ApJ 813:88
############################################################################

def rpitch_phi(a,hp,alpha,eta,r,rp):
    b = (hp/rp)*np.power(rp/r,1.0+eta)*np.power(r,alpha)/np.abs(np.power(r,alpha)-np.power(rp,alpha))

    #b = 1.0/b

#    pitch = d1*np.power(r,-0.04)
    pitch = np.arctan(b)
    return pitch,b

def rpitchspiral_theta(r,a,hp,alpha,eta,rp,x0,y0,xsign=1,ysign=1):
    pitch, b = rpitch_phi(a,hp,alpha,eta,r,rp) 
#    b =0.3    
    r0 = np.sqrt(x0*x0 + y0*y0)
    theta = np.log((r-r0)/a)/b
    theta = np.mod(theta,2.0*np.pi)

    x = xsign*r*np.cos(theta)
    y = ysign*r*np.sin(theta)
    return x,y,theta,pitch,b

def rpitchspiral_thetam(r,m,xsign=1,ysign=1):
    a = m[0]
    hp = m[1]
    alpha = m[2]
    eta = m[3]
    rp = m[4]
    x0 = m[5]
    y0 = m[6]
    
    return rpitchspiral_theta(r,a,hp,alpha,eta,rp,x0,xsign=xsign,ysign=ysign)

def rpitchspiral_x(t,a,hp,alpha,eta,r,rp,x0,xsign=1):
	# Compute b parameter
        pitch,b = rpitch_phi(a,hp,alpha,eta,r,rp)

	# Now compute spiral x position
	return xsign*a*np.exp(b*t)*np.cos(t) + x0

def rpitchspiral_xm(t,m,r,xsign=1):
    a = m[0]
    hp = m[1]
    alpha = m[2]
    eta = m[3]
    rp = m[4]
    x0 = m[5]

    
    return rpitchspiral_x(t,a,hp,alpha,eta,r,rp,x0,ysign=1)

def rpitchspiral_y(t,a,hp,alpha,eta,r,rp,y0,ysign=1):
    # Compute b parameter

    pitch,b = rpitch_phi(a,hp,alpha,eta,r,rp)

    # Now compute spiral y position
    return ysign*a*np.exp(b*t)*np.sin(t) + y0

def rpitchspiral_ym(t,m,r,xsign=1):
    a = m[0]
    hp = m[1]
    alpha = m[2]
    eta = m[3]
    rp = m[4]
    y0 = m[6]

    return rpitchspiral_y(t,a,hp,alpha,eta,r,rp,y0,ysign=ysign)



def find_minimum_t_rpitchspiral(xi,yi, a,hp,alpha,eta,rp,x0,y0, npoints,xsign=1.0,ysign=1.0):
    '''Find the minimum t value for points (xi,yi) given spiral parameters (a,b,x0,y0)'''

    t = np.linspace(0.0,10.0,num=npoints)

    tmin = -1.0
    sepmin = 1.0e30

    for i in range(npoints):

      
	r = np.sqrt(xi*xi + yi*yi)
        x = rpitchspiral_x(t[i], a, hp,alpha,eta,r,rp,x0, xsign=xsign)
        y = rpitchspiral_y(t[i], a, hp,alpha,eta,r,rp,y0, ysign=ysign)

        sep = separation(xi, yi, x, y)

        if(sep<sepmin):
            sepmin = sep
            tmin = t[i]

    return tmin,sepmin


#
# Generate x,y points for an rpitch  spiral curve
#
def generate_rpitchspiral_curve(xbegin,ybegin,xend,yend,a,hp,alpha,eta,rp,x0,y0,xsign=1,ysign=1,nplot=100):
    '''Generate x,y, points of an rpitch spiral'''

    nfind = 1000

    rmin = np.sqrt(xbegin*xbegin + ybegin*ybegin)
    rmax = np.sqrt(xend*xend + yend*yend)
    print rmin, rmax
    r = np.linspace(rmin,rmax,num=nplot)
    xspiral = np.zeros(nplot)
    yspiral = np.zeros(nplot)

    for i in range(nplot):

        xspiral[i], yspiral[i] = rpitchspiral_theta(a,hp,alpha,eta,rp,x0,xsign=xsign,ysign=ysign)
    
    return xspiral,yspiral

def generate_rpitchspiral_curvem(xbegin,xend,ybegin,yend,m,xsign=1,ysign=1,nplot=100):
    '''Generate x,y, points of an rpitch spiral'''
    a = m[0]
    hp = m[1]
    alpha = m[2]
    eta = m[3]
    rp = m[4]
    x0 = m[5]
    y0 = m[6]
  
    return generate_rpitchspiral_curve(xbegin,xend,ybegin,yend,a,hp,alpha,eta,rp,x0,y0,xsign,ysign,nplot)



def get_chisquared_rpitchspiral(x,y,a,hp,alpha,eta,rp,x0,y0,npoints,xsign=1.0,ysign=1.0,sigma=1.0):

    '''Returns the chi-squared of a r-dependent pitch spiral model given arrays x,y
    Assumes uniform errors'''

    tmin = np.zeros(len(x))
    sepmin = np.zeros(len(x))

    for i in range(len(x)):        
        tmin[i], sepmin[i] = find_minimum_t_rpitchspiral(x[i], y[i], a, hp,alpha,eta,rp, x0, y0, npoints,xsign=xsign,ysign=ysign)

    return np.sum(sepmin)/(2.0*len(x)*sigma*sigma)

    
def opt_chisquared_rpitchspiral(m,x,y,npoints,xsign=1.0,ysign=1.0,sigma=1.0,verbose=True):

    a = m[0]
    hp = m[1]
    alpha = m[2]
    eta = m[3]
    rp = m[4]
    x0 = m[5]
    y0 = m[6]

    chimin = get_chisquared_rpitchspiral(x,y,a,hp,alpha,eta,rp,x0,y0,npoints,xsign,ysign,sigma)

    if(m[0]<0.0 or m[1]<0.0 or m[2]<0.0 or m[3]<0.0 or m[4]<0.0): chimin = 1.0e30

    if(verbose):print 'Chi^2 - %4.2f'%chimin,': Parameters  ',' '.join(['%4.2e']*len(m))%tuple(m)

    return chimin




