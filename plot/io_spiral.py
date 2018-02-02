import numpy as np

#
# Functions for I/O etc
#


spiralchoices = ['logarithmic','hyperbolic','power','rpitch']
spiraltexts = ['Logarithmic Spiral', 'Hyperbolic Spiral','Power Spiral' 'r-dependent pitch spiral']
nspiralchoices = len(spiralchoices)
nspiralparams = [4,3,4,8]

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
    
#
# Generate x,y points for a logarithmic spiral curve
#
def generate_logspiral_curve(tmin,tmax,a,b,x0,y0,xsign=1,ysign=1,npoints=100):
    '''Generate x,y, points of a logarithmic spiral'''

    t = np.linspace(tmin,tmax,num=npoints)
    xspiral = np.zeros(npoints)
    yspiral = np.zeros(npoints)

    for i in range(npoints):
        xspiral[i] = logspiral_x(t[i],a,b,x0,xsign=xsign)
        yspiral[i] = logspiral_y(t[i],a,b,y0,ysign=ysign)
    
    return xspiral,yspiral

def generate_logspiral_curvem(tmin,tmax,m,xsign=1,ysign=1,npoints=100):
    '''Generate x,y, points of a logarithmic spiral'''
    m[0] = a
    m[1] = b
    m[2] = x0
    m[3] = y0
    return generate_logspiral_curve(tmin,tmax,a,b,x0,y0,xsign=xsign,ysign=ysign,npoints=npoints)


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
        
        if(sep<sepmin):
            sepmin = sep
            tmin = t[i]
        

    return tmin,sepmin

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
    if(verbose):print chisquared, m
    return chisquared

#
# End of functions for logarithmic spirals
#


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
def generate_hypspiral_curve(tmin,tmax,c,x0,y0,xsign=1,ysign=1,npoints=100):
    '''Generate x,y, points of a hyperbolic spiral'''

    t = np.linspace(tmin,tmax,num=npoints)
    xspiral = np.zeros(npoints)
    yspiral = np.zeros(npoints)

    for i in range(npoints):
        xspiral[i] = hypspiral_x(t[i],c,x0,xsign=xsign)
        yspiral[i] = hypspiral_y(t[i],c,y0,ysign=ysign)
    
    return xspiral,yspiral

def generate_hypspiral_curvem(tmin,tmax,m,xsign=1,ysign=1,npoints=100):
    '''Generate x,y, points of a hyperbolic spiral'''
    m[0] = c
    m[1] = x0
    m[2] = y0

    return generate_logspiral_curve(tmin,tmax,a,b,x0,y0,xsign=xsign,ysign=ysign,npoints=npoints)


# Find the minimum t value for (xi,yi) given spiral parameters (a,b,x0,y0)
# Multiple local minima possible, must be careful

def find_minimum_t_hypspiral(xi,yi, c,x0,y0, npoints,xsign=1.0,ysign=1.0):
    '''Find the minimum t value for points (xi,yi) given spiral parameters (a,b,x0,y0)'''
    
    t = np.linspace(0.0,100.0,num=npoints)
    
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
    if(verbose):print chisquared, m
    return chisquared

#
# End of hyperbolic spiral functions
#


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
    m[0] = a
    m[1] = n
    m[2] = x0

    return powspiral_x(t,a,b,x0,xsign=xsign)

# Same for y
def powspiral_y(t,a,n,y0,ysign=1):    
    return ysign*a*pow(t,n)*np.sin(t) + y0

def powspiral_ym(t,m,ysign=1):
    m[0] = a
    m[1] = n
    m[3] = y0

    return powspiral_y(t,a,b,y0,ysign=ysign)
    
#
# Generate x,y points for a logarithmic spiral curve
#
def generate_powspiral_curve(tmin,tmax,a,n,x0,y0,xsign=1,ysign=1,npoints=100):
    '''Generate x,y, points of a logarithmic spiral'''

    t = np.linspace(tmin,tmax,num=npoints)
    xspiral = np.zeros(npoints)
    yspiral = np.zeros(npoints)

    for i in range(npoints):
        xspiral[i] = logspiral_x(t[i],a,n,x0,xsign=xsign)
        yspiral[i] = logspiral_y(t[i],a,n,y0,ysign=ysign)
    
    return xspiral,yspiral

def generate_powspiral_curvem(tmin,tmax,m,xsign=1,ysign=1,npoints=100):
    '''Generate x,y, points of a logarithmic spiral'''
    m[0] = a
    m[1] = n
    m[2] = x0
    m[3] = y0
    return generate_powspiral_curve(tmin,tmax,a,n,x0,y0,xsign=xsign,ysign=ysign,npoints=npoints)


# Find the minimum t value for (xi,yi) given spiral parameters (a,b,x0,y0)
# Multiple local minima possible, must be careful

def find_minimum_t_powspiral(xi,yi, a,n,x0,y0, npoints,xsign=1.0,ysign=1.0):
    '''Find the minimum t value for points (xi,yi) given spiral parameters (a,b,x0,y0)'''
    
    t = np.linspace(0.0,10.0,num=npoints)
    
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

    chisquared = get_chisquared_powspiral(x,y,a,n,x0,y0,npoints,xsign,ysign,sigma)
    if(verbose):print chisquared, m
    return chisquared

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
    m[0] = a
    m[1] = hp
    m[2] = alpha
    m[3] = eta
    m[4] = rp
    m[5] = x0
    m[6] = y0
    
    return rpitchspiral_theta(r,a,hp,alpha,eta,rp,x0,xsign=xsign,ysign=ysign)

def rpitchspiral_x(t,a,hp,alpha,eta,r,rp,x0,xsign=1):
	# Compute b parameter
        pitch,b = rpitch_phi(a,hp,alpha,eta,r,rp)

	# Now compute spiral x position
	return xsign*a*np.exp(b*t)*np.cos(t) + x0

def rpitchspiral_xm(t,m,r,xsign=1):
    m[0] = a
    m[1] = hp
    m[2] = alpha
    m[3] = eta
    m[4] = rp
    m[5] = x0
    
    return rpitchspiral_x(t,a,hp,alpha,eta,r,rp,x0,ysign=1)

def rpitchspiral_y(t,a,hp,alpha,eta,r,rp,y0,ysign=1):
    # Compute b parameter

    pitch,b = rpitch_phi(a,hp,alpha,eta,r,rp)

    # Now compute spiral y position
    return ysign*a*np.exp(b*t)*np.sin(t) + y0

def rpitchspiral_ym(t,m,r,xsign=1):
    m[0] = a
    m[1] = hp
    m[2] = alpha
    m[3] = eta
    m[4] = rp
    m[6] = y0
    
    return rpitchspiral_x(t,a,hp,alpha,eta,r,rp,y0,ysign=ysign)



    



def find_minimum_t_rpitchspiral(xi,yi, a,hp,alpha,eta,rp,x0,y0, npoints,xsign=1.0,ysign=1.0):
    '''Find the minimum t value for points (xi,yi) given spiral parameters (a,b,x0,y0)'''

    t = np.linspace(0.0,10.0,num=npoints)

    tmin = -1.0
    sepmin = 1.0e30

    for i in range(npoints):

	r = np.sqrt(xi*xi + yi*yi)
        x = rpitchspiral_x(t[i], a, hp,alpha,eta,r,rp,x0, xsign=xsign)
        y = rpitchspiral_y(t[i], a, hp,alpha,eta,r,rp, y0, ysign=ysign)

        sep = separation(xi, yi, x, y)

        if(sep<sepmin):
            sepmin = sep
            tmin = t[i]

    return tmin,sepmin




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

    if(verbose):print chimin, m
    return chimin




