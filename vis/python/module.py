import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 20}

mpl.rc('font', **font) 
mycmap    = 'RdBu_r'

# useful physical constants in cgs unit
c   = 2.99792458e10    #speed of light cm/s
G   = 6.6743e-8        #gravitational constant cm^3g^-1s^-2
h   = 6.6261e-27       #planck constnat cm^2gs^-1
k   = 1.3807e-16       #boltzman constant cm2 g s-2 K-1
sig = 5.6704e-5        #stepen-boltzman constant g s-3 K-4
e   = 4.8032e-10       #electron charge cm3/2 g1/2 s-1
me  = 9.1094e-28       #electron mass (g)
mp  = 1.6726e-24       #proton mass (g)
eV  = 1.602177e-12     #electron volt (erg)

# Astronomical constants in cgs unit
Msun = 1.989e33        #solar mass
Mear = 5.974e27        #earth mass
Mjup = 1.899e30        #jupyter mass
Rsun = 6.955e10        #solar radius
Rear = 6.378e8         #earth radius
Rjup = 7.149e9         #jupyter radius
Lsun = 3.839e33        #solar luminosity (erg/s)
AU   = 1.496e13        #astronomical unit (cm)
Ly   = 9.461e17        #light year
pc   = 3.086e18        #parsec
year = 3.156e7         #year
