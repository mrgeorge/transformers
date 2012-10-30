#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import numpy as np
import scipy.integrate
from esutil import cosmology

cosmo=cosmology.Cosmo(h=0.72,omega_m=0.258,omega_l=0.742) # WMAP5
H0yrInv=cosmo.H0() / (3.08e19) * (3.15e7) # H0 converted to yr-1

# Fakhouri, Ma, & Boylan-Kolchin 2010 Eq 2 
def mdotHalo(redshift,mhalo,medOption):
    if(medOption):
        AA=25.3
        BB=1.65
        alpha=1.1
    else: #mean is default
        AA=46.1
        BB=1.11
        alpha=1.1

    mdot= AA * (mhalo/10**12)**alpha * (1. + BB*redshift) * np.sqrt(cosmo.omega_m() * (1.+redshift)**3 + cosmo.omega_l())

    return mdot


# Fakhouri, Ma, & Boylan-Kolchin 2010 Eq 2 converted to dM/dz
def dmdzHalo(redshift,mhalo,medOption):
    if(medOption):
        AA=25.3
        BB=1.65
        alpha=1.1
    else: #mean is default
        AA=46.1
        BB=1.11
        alpha=1.1

    dmdz= -AA/H0yrInv * (mhalo/10**12)**alpha * (1. + BB*redshift) / (1.+redshift)

    return dmdz

def haloGrowth(mhalo,zi,zf,medOption=False):
    haloGrowth=scipy.integrate.quad(dmdzHalo,zi,zf,args=(mhalo,medOption))[0]/mhalo + 1.

    return haloGrowth # final/initial halo mass

def ageOfUniverse(redshift):
    scale0=1./(1.+redshift)
    age=scipy.integrate.quad(lambda scale: (1./H0yrInv) * (1./scale) / np.sqrt(cosmo.omega_m() * scale**-3 + cosmo.omega_l()), 0, scale0)

    return age[0]/1.e9 #units are Gyr

def stellarGrowth(zi,zf): # note independent of stellar mass
    AA=26.
    BB=-2.2
    stellarGrowth= np.exp((-AA/1.2) * (ageOfUniverse(zf)**(-1.2) - ageOfUniverse(zi)**(-1.2)))

    return stellarGrowth # final/initial stellar mass

# MAIN - if called from command line
if __name__ == '__main__':

    # how much does a 10**13.5 mass halo grow from z=1 to 0, and z=0.9 to 0.35?

    mh=10.**13.5

    print np.log10(haloGrowth(mh,1,0)) # growth in dex
    print np.log10(haloGrowth(mh,0.9,0.35))


    # how much does a star-forming galaxy grow over the same redshift ranges?
    
    print np.log10(stellarGrowth(1,0))
    print np.log10(stellarGrowth(0.9,0.35))

    print np.log10(stellarGrowth(0.9,0.65))
    print np.log10(stellarGrowth(0.65,0.35))
