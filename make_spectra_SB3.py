# # # # # # # Simple script to create "mock SB1/SB2" data  # # # # # # # 
# Created by Tomer Shenar, T.Shenar@uva.nl; tomer.shenar@gmail.com


import matplotlib.pyplot as plt
import glob
import os
import numpy as np
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

T0_out = 34.
P_out = 200.
e_out = 0.3
omega_out = 160. * np.pi/180.
K_out = 29.

#reflex motion of binary
K_bin = ( P/P_out * ( (1 - e**2)/(1 - e_out**2) )**(3./2.) * (K1 + K2)**3 / K_out )**0.5 - K_out

print("reflex binary K: ", K_bin)

if K_bin<0:
    print('Configuration not possible, K1 + K2 must fulfil:')
    print('K1 + K2 > ', K_out * ( P_out / P * ((1 - e_out**2)/(1 - e**2) )**(3./2.) )**(1./3.) )
    sys.exit()



# Data properties (band, resolution and sampling, i.e. Dlam, & S2N)
lamB = 4000.
lamR = 5000.
Resolution = 20000.
SamplingFac = 3
Dlam = (lamB+lamR)/2. /Resolution / SamplingFac
S2N = 100. 

#Flux ratio, F2_to_ftot F3_to_ftot 
Q2 = 0.3
Qout = 0.3



# Number of epochs
specnum = 40

## Mask path
MaskPath = '/Users/tomer/Desktop/Work/models/TLUSTY/G40000g400v10.vis.rectvmac30vsini200.dat'
MaskPath2 = '/Users/tomer/Desktop/Work/models/TLUSTY/BG20000g300v2.vis.rectvmac30vsini300.dat'
MaskPath3 = '/Users/tomer/Desktop/Work/models/TLUSTY/G35000g400v10.vis.recvmac30vsini100.dat'

#MaskPath = '/Users/tomer/Desktop/Work/models/TLUSTY/BG27000g475v2.vis.rect'
#MaskPath2 = '/Users/tomer/Desktop/Work/models/TLUSTY/BG20000g300v2.vis.rect'
#MaskPath3 = '/Users/tomer/Desktop/Work/models/TLUSTY/G35000g400v10.vis.rect'



# Nebular Lines?
NebularLines = False
Nebwaves = np.array([4026., 4102., 4120.8, 4143., 4340.5,  4388, 4471.5])

######################################
# # # # # # #END USER INPUT # # # # # 
######################################


clight = 2.9979E5

efac = np.sqrt((1 + e) / (1 - e))
efac_out = np.sqrt((1 + e_out) / (1 - e_out))

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
MaskTemp3 = np.loadtxt(MaskPath3)

Waves1 = MaskTemp[:,0] + np.random.normal(0,1E-10, len(MaskTemp))
Waves2 = MaskTemp2[:,0] + np.random.normal(0,1E-10, len(MaskTemp2))
Waves3 = MaskTemp3[:,0] + np.random.normal(0,1E-10, len(MaskTemp3))

Mask = interp1d(Waves1, MaskTemp[:,1], bounds_error=False, fill_value=1., kind='cubic')(wavegrid)
Mask2 = interp1d(Waves2,MaskTemp2[:,1], bounds_error=False, fill_value=1., kind='cubic')(wavegrid)
Mask3 = interp1d(Waves3,MaskTemp3[:,1], bounds_error=False, fill_value=1., kind='cubic')(wavegrid)


stdConv =  DlamRes/Dlam /np.sqrt(2*np.log(2))/2.
kernel = Gaussian1DKernel(stddev=stdConv)

Mask = convolve(Mask, kernel, normalize_kernel=True, boundary='extend') 
Mask2 = convolve(Mask2, kernel, normalize_kernel=True, boundary='extend') 
Mask3 = convolve(Mask3, kernel, normalize_kernel=True, boundary='extend') 

np.savetxt('obs/Atemp.txt', np.c_[wavegrid, Mask])
np.savetxt('obs/Btemp.txt', np.c_[wavegrid, Mask2])
np.savetxt('obs/Ctemp.txt', np.c_[wavegrid, Mask3])


sig = 1/S2N      

#Nebula
if NebularLines:
    Mask3_pre = np.copy(wavegrid) * 0. + 0.
    IndsNebWaves = np.array([np.argmin(np.abs(wavegrid-Nebwave)) for Nebwave in Nebwaves])
    Mask3_pre[IndsNebWaves] += NebStrengh
    MaskNeb = convolve(Mask3_pre, kernel, normalize_kernel=True,
                        boundary='extend') 
else:
    MaskNeb = Mask*0.

phasesfile.write('MJD obsname\n')
for i in np.arange(specnum):
      NebFac = 1.
      MJD = random.uniform(0., 1.)*max(P, P_out)
      phase = (MJD - T0)/P - int((MJD-T0)/P)
      phase_out = (MJD - T0_out)/P_out - int((MJD-T0_out)/P_out)
      M = 2 * np.pi * phase
      M_out = 2 * np.pi * phase_out
      E =  Kepler(1., M, e)
      E_out =  Kepler(1., M_out, e_out)
      nu = 2. * np.arctan(efac * np.tan(0.5 * E))     
      nu_out = 2. * np.arctan(efac_out * np.tan(0.5 * E_out))     
      v1, v2= v1v2(nu, Gamma, K1, K2, omega, e)
      vTertiary, vBin = v1v2(nu_out, Gamma, K_out, K_bin, omega_out, e_out)
      Facshift1 = np.sqrt( (1 + (v1 + vBin)/clight) / (1 - (v1 + vBin)/clight))
      Facshift2 = np.sqrt( (1 + (v2 + vBin)/clight) / (1 - (v2 + vBin)/clight))      
      Facshift1_out = np.sqrt( (1 + vTertiary/clight) / (1 - vTertiary/clight))      
      Maskshift1 = interp1d(wavegrid*Facshift1, Mask, bounds_error=False, fill_value=1., kind='cubic')(wavegrid)
      Maskshift2 = interp1d(wavegrid*Facshift2, Mask2, bounds_error=False, fill_value=1., kind='cubic')(wavegrid)   
      Maskshift3 = interp1d(wavegrid*Facshift1_out, Mask3, bounds_error=False, fill_value=1., kind='cubic')(wavegrid)   
      MaskSums = (1-Q2 - Qout)*Maskshift1 + Q2*Maskshift2 + Qout*Maskshift3 + NebFac*MaskNeb 
      #convoluted = convolve(MaskSums, kernel, normalize_kernel=True, boundary='extend')
      noiseobs = MaskSums + np.random.normal(0,sig, len(wavegrid))
      obsname = 'obs/obs_' + str(i) +  "_V1_" + str(v1) + "_V2_" + str(v2) + "_V3_" + str(vTertiary)
      np.savetxt(obsname, np.c_[wavegrid, noiseobs])
      phasesfile.write(str(MJD) + ' '  + obsname + '\n')
      #if i==0:
          #np.savetxt('Observation_example.txt', np.c_[wavegrid, noiseobs])

phasesfile.close()      
