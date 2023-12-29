# LISA PSD
 


#############################################################################
###################                                       ###################
###################             LISA PSD                  ###################
###################                                       ###################
#############################################################################



import numpy as np
from scipy import interpolate


fm     = 3.168753575e-8   # LISA modulation frequency
YEAR   = 3.15581497632e7  # year in seconds
AU     = 1.49597870660e11 # Astronomical unit (meters)
Clight = 299792458.       # speed of light (m/s)

Tobs=4*YEAR
Larm=2.5e9
NC=2

fstar = Clight/(2*np.pi*Larm)

def Pn(f):
    """
    Calculate the Strain Power Spectral Density
    """
    
    # single-link optical metrology noise (Hz^{-1}), Equation (10)
    P_oms = (1.5e-11)**2*(1. + (2.0e-3/f)**4) 
    
    # single test mass acceleration noise, Equation (11)
    P_acc = (3.0e-15)**2*(1. + (0.4e-3/f)**2)*(1. + (f/(8.0e-3))**4) 
    
    # total noise in Michelson-style LISA data channel, Equation (12)
    Pn = (P_oms + 2.*(1. + np.cos(f/fstar)**2)*P_acc/(2.*np.pi*f)**4)/Larm**2
    
    return Pn
    
def SnC(f):
    """
    Get an estimation of the galactic binary confusion noise are available for
        Tobs = {0.5 yr, 1 yr, 2 yr, 4yr}
    Enter Tobs as a year or fraction of a year
    """

    # Fix the parameters of the confusion noise fit
    if (Tobs < .75*YEAR):
        est = 1
    elif (0.75*YEAR < Tobs and Tobs < 1.5*YEAR):
        est = 2
    elif (1.5*YEAR < Tobs and Tobs < 3.0*YEAR):   
        est = 3
    else:
        est = 4
            
    if (est==1):
        alpha  = 0.133
        beta   = 243.
        kappa  = 482.
        gamma  = 917.
        f_knee = 2.58e-3  
    elif (est==2):
        alpha  = 0.171
        beta   = 292.
        kappa  = 1020.
        gamma  = 1680.
        f_knee = 2.15e-3 
    elif (est==3):
        alpha  = 0.165
        beta   = 299.
        kappa  = 611.
        gamma  = 1340.
        f_knee = 1.73e-3  
    else:
        alpha  = 0.138
        beta   = -221.
        kappa  = 521.
        gamma  = 1680.
        f_knee = 1.13e-3 
    
    A = 1.8e-44/NC
    
    Sc  = 1. + np.tanh(gamma*(f_knee - f))
    Sc *= np.exp(-f**alpha + beta*f*np.sin(kappa*f))
    Sc *= A*f**(-7./3.)
    
    return Sc
    
def Sn(f):
    """ Calculate the sensitivity curve """
    # Proper sky averaged response function (explained a bit here: https://arxiv.org/pdf/1803.01944.pdf)
    R = interpolate.splev(f, R_INTERP, der=0)
    Sn = Pn(f)/R + SnC(f)

    return Sn

transfer_data = np.genfromtxt("R.txt")
        
f = transfer_data[:,0]*fstar        # convert to frequency
R = transfer_data[:,1]*NC           # response gets improved by more data channels
R_INTERP = interpolate.splrep(f, R, s=0)



def PSDfunc(f):
    # PSD for LISA
    highcutoff = 2. # "hard-coded"
    f1 = np.array([i for i in f if i <= highcutoff])
    f2 = np.array([i for i in f if i > highcutoff])
    if len(f2) > 0:
        P = np.concatenate(([Sn(f[1])],Sn(f1[1:]),Sn(f[len(f)-len(f2)])*np.ones(len(f2))))
    else:
        P = np.concatenate(([Sn(f[1])],Sn(f1[1:])))
    return P


