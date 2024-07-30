# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 21:05:37 2018

@author: dselva
"""
import numpy as np
c = 3e8
from orbits import R_earth

def maxPRF(hkm,eta):
    r = 1000*hkm/np.cos(np.radians(eta))
    return c/2/r
    
def range_res(B,theta):
    return c/2/B/np.sin(np.radians(theta))

def azimuth_res(hkm,f,L,theta):
    from comms import f2l
    h = 1000*hkm
    return h*f2l(f)/L/np.cos(np.radians(theta))

def focused_SAR_azimuth_res(L):
    return L/2

def unfocused_SAR_azimuth_res(f,hkm):
    from comms import f2l
    h = 1000*hkm
    return np.sqrt(2*h*f2l(f))

def size_radar(range_res,azimuth_res,hkm,f,theta):
    from comms import f2l
    B = c/2/range_res/np.sin(np.radians(theta))
    h = 1000*hkm
    L = h*f2l(f)/np.cos(np.radians(theta))/azimuth_res
    return L,B

def size_SAR(range_res,azimuth_res,hkm,f,theta):
    B = c/2/range_res/np.sin(np.radians(theta))
    L = 2*azimuth_res
    return L,B

def radar_equation(Pt,L,W,hkm,f,B,T,theta,sigma_dB,hrz_avg=False,SAR=False,apert_type="rect"):
    from comms import f2l
    from remote_sensing import kb
    
    l = f2l(f)
    if apert_type == 'rect':
        eff = 0.55
        Aeff = eff*L*W
    elif apert_type == 'circ':
        eff = 0.55
        Aeff = eff*np.pi*L**2/4
        
    G = 4*np.pi*Aeff/l**2
    Xr= range_res(B,theta)
    
    if SAR == 'Unfocused':
        Xa = unfocused_SAR_azimuth_res(f,hkm)
    elif SAR == 'Focused':
        Xa = focused_SAR_azimuth_res(L)
    else:
        Xa = azimuth_res(hkm,f,L,theta)
    
    if hrz_avg>0:
        Nr = hrz_avg/Xr
        if Nr<1:
            Nr = 1
        Na= hrz_avg/Xa
        if Na<1:
            Na=1
    else:
        Na=1
        Nr=1
    area = Na*Xa*Nr*Xr*np.cos(np.radians(theta))*(10**(sigma_dB/10))
    r = 1000*hkm/np.cos(np.radians(theta))
    S=(Pt*G/4/np.pi/r**2)*area*(Aeff/4/np.pi/r**2)
    SNR = S/(kb*T*B)
    return SNR

def SAR_equation(Pt,tau,PRF,L,W,hkm,f,B,T,theta,sigma_dB,SAR="Focused",apert_type="rect",RE=R_earth):
    from comms import f2l
    from remote_sensing import kb
    from orbits import orb_vel_c
    l = f2l(f)
    if apert_type == 'rect':
        eff = 0.55
        Aeff = eff*L*W
    elif apert_type == 'circ':
        eff = 0.55
        Aeff = eff*np.pi*L**2/4
        
    G = 4*np.pi*Aeff/l**2
    Xr= range_res(B,theta)
    
    if SAR == 'Unfocused':
        Xa = unfocused_SAR_azimuth_res(f,hkm)
    elif SAR == 'Focused':
        Xa = focused_SAR_azimuth_res(L)
    else:
        Xa = azimuth_res(hkm,f,L,theta)
     
    Gr = tau*B
    v = RE*orb_vel_c(hkm)/(RE+hkm)
    Ga = PRF*l*hkm*1000/Xa/v
    
    area = Gr*Ga*Xa*Xr*np.cos(np.radians(theta))*(10**(sigma_dB/10))
    r = 1000*hkm/np.cos(np.radians(theta))
    S=(Pt*G/4/np.pi/r**2)*area*(Aeff/4/np.pi/r**2)
    SNR = S/(kb*T*B)
    return SNR


def noise_equiv_sigma(Pt,L,W,hkm,f,B,T,theta,sigma_dB,hrz_avg=False,SAR=False):
    SNR= radar_equation(Pt,L,W,hkm,f,B,T,theta,0,hrz_avg,SAR)
    sigma = 1/SNR
    return 10*np.log10(sigma)
    
def size_radar2(Pt,L,hkm,f,B,T,theta,sigma_dB,hrz_avg=False,SAR=False):

    from scipy.optimize import minimize
    err = lambda W:(abs(100*np.log10(radar_equation(Pt,L,W,hkm,f,B,T,theta,sigma_dB,hrz_avg,SAR)) - 0))
    W0=1
    res = minimize(err,W0)
    zeW = res.x[0]        
    #Xa = focused_SAR_azimuth_res(L)
    if SAR == 'Unfocused':
        Xa = unfocused_SAR_azimuth_res(f,hkm)
    elif SAR == 'Focused':
        Xa = focused_SAR_azimuth_res(L)
    else:
        Xa = azimuth_res(hkm,f,L,theta)
        
    Xr = range_res(B,theta)
    
    if hrz_avg>0:
        Nr = hrz_avg/Xr
        if Nr<1:
            Nr = 1
        Na= hrz_avg/Xa
        if Na<1:
            Na=1
    else:
        Na=1
        Nr=1
        
    Xa = Na*Xa
    Xr = Nr*Xr    
    return zeW,L,Xa,Xr