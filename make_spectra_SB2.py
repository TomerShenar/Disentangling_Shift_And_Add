# Short script to create mock SB2 binaries, by Tomer Shenar


import matplotlib.pyplot as plt
import glob
import os
import numpy as np
# # # # # # # Simple script to create "mock SB1/SB2" data  # # # # # # # 
# Created by Tomer Shenar, T.Shenar@uva.nl; tomer.shenar@gmail.com

import sys
from astropy.io import ascii
from scipy.interpolate import interp1d
from scipy import stats  
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
import random
from astropy.convolution import Gaussian1DKernel
from astropy.convolution import convolve

######################################
# # # # # # # USER INPUT # # # # # # #
######################################



#Orbit Pars:
T0 = 0.
P = 18.
e = 0.3
omega = 60. * np.pi/180.
K1 = 87.
K2 = 135.
Gamma = 0.


# Data properties (band, resolution and sampling, i.e. Dlam, & S2N)
lamB = 4000.
lamR = 5000.
Resolution = 15000.
SamplingFac = 3
Dlam = (lamB+lamR)/2. /Resolution / SamplingFac
S2N = 100. 

#Flux ratio, F2_to_ftot 
Q = 0.3

# Number of epochs
specnum = 20

# Mask path
MaskPath = '/Users/tomer/Desktop/Work/models/TLUSTY/G40000g400v10.vis.rectvmac30vsini200.dat'
MaskPath2 = '/Users/tomer/Desktop/Work/models/TLUSTY/BG27000g475v2.vis.rectvmac30vsini25.dat'



# Nebular Lines?
NebularLines = False
Nebwaves = np.array([4026., 4102., 4120.8, 4143., 4340.5,  4388, 4471.5])

######################################
# # # # # # #END USER INPUT # # # # # 
######################################


clight = 2.9979E5

efac = np.sqrt((1 + e) / (1 - e))

def v1(nu, Gamma, K1, omega, ecc):  
      v1 = Gamma + K1*(np.cos(omega + nu) + ecc* np.cos(omega))
      #v2 = Gamma + K2*(np.cos(np.pi + Omega + nu) + ecc* np.cos(np.pi + Omega))
      #return np.column_stack((v1,v2))
      return v1

# For converting Mean anomalies to eccentric anomalies (M-->E)
def Kepler(E, M, ecc):
      E2 = (M - ecc*(E*np.cos(E) - np.sin(E))) / (1. - ecc*np.cos(E))
      eps = np.abs(E2 - E) 
      if np.all(eps < 1E-10):
            return E2
      else:
            return Kepler(E2, M, ecc)
        
def v1v2(nu, Gamma, K1, K2, omega, ecc):  
      v1 = Gamma + K1*(np.cos(omega + nu) + ecc* np.cos(omega))
      v2 = Gamma + K2*(np.cos(np.pi + omega + nu) + ecc* np.cos(np.pi + omega))
      return v1, v2
        

#files = glob.glob('/YOUR/PATH/*')
for f in glob.glob('obs/obs*'):
    os.remove(f)
for f in glob.glob('obs/ObsDat.t*'):
    os.remove(f)

### Path to observations: Directory in which the spectra are found:
PathToObservations = '../'


### Path to output: Directory in which output should be written:
PathToOutput = "./"



if not os.path.exists('obs/'):
    os.mkdir('obs/')


phasesfile = open('obs/ObsDat.txt', 'w')

#lamB, lamR = np.loadtxt(MaskPath)[:,0][[0,-1]]
lammid = (lamB + lamR)/2.
DlamRes = lammid/Resolution


wavegrid = np.arange(lamB, lamR, Dlam)

stdConv =  DlamRes/Dlam /np.sqrt(2*np.log(2))/2.
kernel = Gaussian1DKernel(stddev=stdConv)


MaskTemp = np.loadtxt(MaskPath)
MaskTemp2 = np.loadtxt(MaskPath2)

Waves1 = MaskTemp[:,0] + np.random.normal(0,1E-10, len(MaskTemp))
Waves2 = MaskTemp2[:,0] + np.random.normal(0,1E-10, len(MaskTemp2))

Mask = interp1d(Waves1, MaskTemp[:,1], bounds_error=False, fill_value=1., kind='cubic')(wavegrid)
Mask2 = interp1d(Waves2,MaskTemp2[:,1], bounds_error=False, fill_value=1., kind='cubic')(wavegrid)


stdConv =  DlamRes/Dlam /np.sqrt(2*np.log(2))/2.
kernel = Gaussian1DKernel(stddev=stdConv)

Mask = convolve(Mask, kernel, normalize_kernel=True, boundary='extend') 
Mask2 = convolve(Mask2, kernel, normalize_kernel=True, boundary='extend') 

np.savetxt('obs/Atemp.txt', np.c_[wavegrid, Mask])
np.savetxt('obs/Btemp.txt', np.c_[wavegrid, Mask2])


sig = 1/S2N      

#Nebula
if NebularLines:
    Mask3_pre = np.copy(wavegrid) * 0. + 0.
    IndsNebWaves = np.array([np.argmin(np.abs(wavegrid-Nebwave)) for Nebwave in Nebwaves])
    Mask3_pre[IndsNebWaves] += NebStrengh
    Mask3 = convolve(Mask3_pre, kernel, normalize_kernel=True,
                        boundary='extend') 
else:
    Mask3 = Mask*0.

phasesfile.write('MJD obsname\n')
for i in np.arange(specnum):
      NebFac = 1.
      MJD = random.uniform(0., 1.)*P
      phase = (MJD - T0)/P - int((MJD-T0)/P)
      M = 2 * np.pi * phase
      E =  Kepler(1., M, e)
      nu = 2. * np.arctan(efac * np.tan(0.5 * E))     
      v1, v2= v1v2(nu, Gamma, K1, K2, omega, e)
      Facshift1 = np.sqrt( (1 + v1/clight) / (1 - v1/clight))
      Facshift2 = np.sqrt( (1 + v2/clight) / (1 - v2/clight))      
      Maskshift1 = interp1d(wavegrid*Facshift1, Mask, bounds_error=False, fill_value=1., kind='cubic')(wavegrid)
      Maskshift2 = interp1d(wavegrid*Facshift2, Mask2, bounds_error=False, fill_value=1., kind='cubic')(wavegrid)   
      MaskSums = (1-Q)*Maskshift1 + Q*Maskshift2 + NebFac*Mask3 
      #convoluted = convolve(MaskSums, kernel, normalize_kernel=True, boundary='extend')
      noiseobs = MaskSums + np.random.normal(0,sig, len(wavegrid))
      obsname = 'obs/obs_' + str(i) +  "_V1_" + str(v1) + "_V2_" + str(v2)
      np.savetxt(obsname, np.c_[wavegrid, noiseobs])
      phasesfile.write(str(MJD) + ' '  + obsname + '\n')
      #if i==0:
          #np.savetxt('Observation_example.txt', np.c_[wavegrid, noiseobs])

phasesfile.close()      
