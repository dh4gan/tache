import numpy as np

# Functions to deliver parametric forms of spiral (x) (and derivatives)
def logspiral_x(t,a,b,x0,xsign=1):    
    return xsign*a*np.exp(b*t)*np.cos(t) + x0

def logspiral_dx(t,a,b,x0,xsign=1):    
    return -xsign*b*a*np.exp(b*t)*np.sin(t)

def logspiral_d2x(t,a,b,x0,xsign=1):    
    return -xsign*b*b*a*np.exp(b*t)*np.cos(t)

# Same for y

def logspiral_y(t,a,b,y0,ysign=1):    
    return ysign*a*np.exp(b*t)*np.sin(t) + y0

def logspiral_dy(t,a,b,y0,ysign=1):    
    return ysign*b*a*np.exp(b*t)*np.cos(t)

def logspiral_d2y(t,a,b,y0,ysign=1):    
    return -ysign*b*b*a*np.exp(b*t)*np.sin(t)


# Functions for spiral with r-dependent pitch angle

def rpitch_phi(a,d1,d2,eta,r,rp):
    b = d1*np.power(rp/r,1.0+eta)*np.power(r,d2)/np.abs(np.power(r,d2)-np.power(rp,d2))

#    pitch = d1*np.power(r,-0.04)
    pitch = np.arctan(b)
    return pitch,b

def rpitchspiral_theta(r,a,d1,d2,eta,rp,x0,y0,xsign=1,ysign=1):
    pitch, b = rpitch_phi(a,d1,d2,eta,r,rp) 
#    b =0.3    
    r0 = np.sqrt(x0*x0 + y0*y0)
    theta = np.log(r-r0/a)/b
#    theta = np.mod(theta,2.0*np.pi)

    x = xsign*r*np.cos(theta)
    y = ysign*r*np.sin(theta)
    return x,y,theta,pitch,b


def rpitchspiral_x(t,a,d1,d2,eta,r,rp,x0,xsign=1):
	# Compute b parameter
        pitch,b = rpitch_phi(a,d1,d2,eta,r,rp)

	# Now compute spiral x position
	return xsign*a*np.exp(b*t)*np.cos(t) + x0

def rpitchspiral_y(t,a,d1,d2,eta,r,rp,y0,ysign=1):
        # Compute b parameter

	pitch,b = rpitch_phi(a,d1,d2,eta,r,rp)

        # Now compute spiral y position
        return ysign*a*np.exp(b*t)*np.sin(t) + y0


def separation(x1,y1,x2,y2):
    '''Return separation of x y coordinates'''
    return np.sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))
    


# Find the minimum t value for points (xi,yi) given spiral parameters (a,b,x0,y0)
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

def find_minimum_t_rpitchspiral(xi,yi, a,d1,d2,eta,rp,x0,y0, npoints,xsign=1.0,ysign=1.0):
    '''Find the minimum t value for points (xi,yi) given spiral parameters (a,b,x0,y0)'''

    t = np.linspace(0.0,10.0,num=npoints)

    tmin = -1.0
    sepmin = 1.0e30

    for i in range(npoints):

	r = np.sqrt(xi*xi + yi*yi)
        x = rpitchspiral_x(t[i], a, d1,d2,eta,r,rp,x0, xsign=xsign)
        y = rpitchspiral_y(t[i], a, d1,d2,eta,r,rp, y0, ysign=ysign)

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


def opt_chisquared_logspiral(m,x,y,npoints,xsign=1.0,ysign=1.0,sigma=1.0):
    '''Wrapper for scipy.optimize: the chi-squared of a logarithmic spiral model given arrays x,y
    Assumes uniform errors'''
    
    a = m[0]
    b = m[1]
    x0 = m[2]
    y0 = m[3]
  
    chisquared = get_chisquared_logspiral(x,y,a,b,x0,y0,npoints,xsign,ysign,sigma)
    print chisquared, m
    return chisquared

def get_chisquared_rpitchspiral(x,y,a,d1,d2,eta,rp,x0,y0,npoints,xsign=1.0,ysign=1.0,sigma=1.0):

    '''Returns the chi-squared of a r-dependent pitch spiral model given arrays x,y
    Assumes uniform errors'''

    tmin = np.zeros(len(x))
    sepmin = np.zeros(len(x))

    for i in range(len(x)):        
        tmin[i], sepmin[i] = find_minimum_t_rpitchspiral(x[i], y[i], a, d1,d2,eta,rp, x0, y0, npoints,xsign=xsign,ysign=ysign)

    return np.sum(sepmin)/(2.0*len(x)*sigma*sigma)

    
def opt_chisquared_rpitchspiral(m,x,y,npoints,xsign=1.0,ysign=1.0,sigma=1.0):

    a=m[0]
    d1 = m[0]
    d2 = m[1]
    eta = m[2]
    rp = m[3]
    x0 = m[4]
    y0 = m[5]

    chimin = get_chisquared_rpitchspiral(x,y,a,d1,d2,eta,rp,x0,y0,npoints,xsign,ysign,sigma)

    print chimin, m
    return chimin
