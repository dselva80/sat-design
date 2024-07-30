# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 21:03:39 2018

@author: Dani Selva
"""

import matplotlib.pyplot as plt
import numpy as np

h = 6.626e-34
kb = 1.38e-23
c = 3e8

### Planck
def planck(lamda,T,draw=False):
    "Returns the spectral radiance at wavelength lamda for a blackbody at temperature T in W/m^2/Sr/\mum"
    b = 1e-6*2*h*c**2/(lamda**5) * 1 / (np.exp((h*c)/(lamda*kb*T))  - 1)
    ## plot
    if draw:
        plt.plot(1e6*lamda,b)
        plt.xlabel("wavelength (\mum)")
        plt.ylabel("spectral radiance (W/m^2/Sr/\mum)")
    return b

def inv_planck(lamda,b):
    """Returns the temperature in K corresponding to a given radiance in W/m^2/Sr/\mum at wavelength lamda"""
    from scipy.optimize import minimize
    
    err = lambda T:(np.abs(planck(lamda,T) - b))
    T0=300
    res = minimize(err,T0)
    return res.x

def planck_f(f,T,draw=False):
    """Returns spectral radiance in W/m^2 HzSr"""
    b = 2*h*f**3/c**2 * (np.exp(h*f/kb/T) - 1)**(-1)
    if draw:
        plt.plot(1e-9*f,b)
        plt.xlabel("frequency (GHz)")
        plt.ylabel("spectral radiance (W/m^2/Sr/Hz)")
    return b

### Gray bodies
def brightnessT(Tphys,e,f):
    Tb=(kb/h/f*np.log(1+(np.exp(h*f/kb/Tphys)-1)/e))**(-1)
    return Tb
      
### Stefan Boltzmann, Wien's laws
def stefan_boltzmann(e,T,A=1,Tcold=0):
    sigma = 5.67e-8
    Q = sigma*e*(T**4 - Tcold**4)*A
    return Q

def wien(T):
    """Returns the wavelength in um of maximum radiation for a blackbody at temperature T"""
    return 2898/T     
    
### Rayleigh Jeans
def rayleigh_B(f,T):
    return 2*kb*T*f**2/(c**2)
    
def rayleigh_TB(e,T):
    return e*T
    
## Received quantities
def scale_irradiance(I1,R1,r):
    return I1*(R1/r)**2

def received_power_IR(l1,l2,T,pixel_area,r,apert):
   delta_lamda = l2 - l1 # in meters
   l0 = (l1 + l2)/2
   L = planck(l0,T) # spectral radiance in W/m^2 um Sr
   omega = pixel_area / r**2 # solid angle subtended by the emitting pixel from the satellite, in Sr
   Pr = L*apert*omega*delta_lamda*1e6 # radiant flux in W, assuming delta_lamda is small ,and pixel is small compared to distance and radiates isotropically 
   return Pr

def Pr2T(P,l1,l2,pixel_size,rkm,apert):
    delta_lamda = 1e6*(l2 - l1) # in meters
    l0 = (l1 + l2)/2
    L = P/(delta_lamda*np.pi*(apert**2)/4*(pixel_size**2)/(1000*rkm)**2)
    print (delta_lamda,l0,L) 
    T = inv_planck(l0,L)
    return T