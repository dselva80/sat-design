# -*- coding: utf-8 -*-
"""
Created on Fri May 11 14:35:37 2018

@author: Dani Selva
"""
import numpy as np
from orbits import R_earth

### Utilities
def dB(lin):
    return 10*np.log10(lin)

def lin(dB):
    return 10**(dB/10)

def f2l(f):
    return 3e8/f

def l2f(l):
    return 3e8/l

### Gain
def D2G(D,f,eta):
    lamda = f2l(f)
    g = eta*(np.pi*D/lamda)**2
    return g,dB(g)

def G2D(GdB,f,eta):
    lamda = f2l(f)
    G = lin(GdB)
    D = lamda/np.pi*np.sqrt(G/eta)
    return D

def Aeff_circ(eta,D):
    return eta*np.pi*D**2/4
## Losses

def fsl(f,hkm,elev,RE=R_earth):
    from geometry import comm_range
    r = 1000*comm_range(hkm,elev,RE)
    Ls = 20*np.log10(f2l(f)/(4*np.pi*r))
    return Ls,lin(Ls)
    
def atm_losses(f):
    fGHz = f/1e9
    # Gasses
    if fGHz < 10:
        Latm = 0
    elif fGHz < 30:
        Latm = -1
    elif fGHz < 40:
        Latm = -5
    elif fGHz < 50:
        Latm = -10
    else:
        Latm = -100

    
    L = Latm 
    return L,lin(L)

def get_k_alpha(f):
    f = f/1e9
    if f <2:
        k = 0.0000387
        alpha = 0.912
    elif f < 4:
        k = 0.000154
        alpha = 0.963
    elif f < 10:
        k = 0.00454
        alpha = 1.327
    elif f < 12:
        k = 0.0101
        alpha = 1.276
    elif f < 15:
        k = 0.0188
        alpha = 1.217
    elif f<25:
        k = 0.0751
        alpha = 1.099
    elif f < 30:
        k = 0.124
        alpha = 1.061
    elif f < 35:
        k = 0.187
        alpha = 1.021
    elif f < 50:
        k = 0.35
        alpha = 0.939
    else:
        exit("I don't have data for this f")
    return k,alpha
        
def rain_losses(f,elev,lat=40,RE=R_earth):
    if lat>23:
        h_rain = 5 - 0.075*(lat - 23)
    elif lat<23 and lat>0:
        h_rain = 5
    elif lat > -21 and lat < 0:
        h_rain = 5
    elif lat>-71:
        h_rain = 5 + 0.1*(lat + 21)
    D_rain = h_rain/np.sin(np.radians(elev))
    
    k,alpha = get_k_alpha(f)
    R = 6
    gamma_rain = k*R**alpha
    Lrain = -D_rain*gamma_rain
    return Lrain,lin(Lrain)
    
def antenna_noise_temp(f):
    """Values for DL"""
    import sys
    fGHz = f/1e9
    if fGHz < 1:
        Ta = 150
    elif fGHz < 12:
        Ta = 25
    elif fGHz <= 30:
        Ta = 100
    else:
        sys.exit("Don't have data to estimate antenna noise temperature")
    return Ta

def system_noise_temp(Ta,FdB,LldB=0):
    Ll = lin(LldB)
    Ts = Ta + 290*(lin(FdB)-1)/Ll + 290*(1-Ll)/Ll
    return Ts

def G_T2T(G_TdB,D,f):
    G_T = lin(G_TdB)
    G = D2G(D,f,0.55)[0]
    T= G/G_T
    return T

### Link budget
def EbN0(Pt,Gt,Gr,f,hkm,elev,Ta,FdB,Rb,LldB=0,RE=R_earth,lat = 40):
    
    k = 1.38e-23
    Ls=fsl(f,hkm,elev,RE)[1]
    if Ta == -1:
        Ta=antenna_noise_temp(f)
    Ts = system_noise_temp(Ta,FdB,LldB)
    Latm = atm_losses(f)[1]
    Lrain = rain_losses(f,elev,lat,RE)[1]
    EbN0 = Pt*Gt*Gr*Ls*Latm*Lrain*lin(LldB)/(k*Ts*Rb)
    #print ('EbN0 ',EbN0)
    print(('Pt = {}, Gt(lin) = {}, Gr(lin) = {}, Ls(dB) = {}, Latm(dB) = {}, Lrain(dB) = {}, Ll(dB) = {}, Ts = {}, Rb = {}, EbN0(dB) = {}').format(Pt,Gt,Gr,dB(Ls),dB(Latm),dB(Lrain),LldB,Ts,Rb,dB(EbN0)))
    return EbN0,dB(EbN0)

def min_EbN0_mod(BER,mod):
    from scipy.stats import norm
    if mod == 'QPSK' or mod == 'BPSK':
        EbN0 = (norm.isf(BER)**2)/2
    elif mod == '8PSK':
        EbN0 = ((norm.isf(0.5*3*BER)/np.sin(np.pi/8))**2)/(2*3)
    elif mod == '16QAM':
        EbN0 = 1.25*((norm.isf((4/4)*BER))**2)
    elif mod == '64QAM':
        EbN0 = (63/18)*((norm.isf((6/4)*BER))**2)
    else:
        print("Modulation not returned")
    return EbN0,dB(EbN0)

def min_EbN0_shannon(Rb,B):
    eta = Rb/B
    Eb_N0 = (2**eta - 1)/eta
    return Eb_N0,dB(Eb_N0)
    
def link_margin(Pt,Dt,Grlin,f,hkm,elev,Ta,FdB,Rb,B,BER,mod = 'QPSK',LldB=0,RE=R_earth):
    print ('Dt ', Dt)
    Gt = D2G(Dt,f,0.55)[0]
    minEbN0 = max(min_EbN0_shannon(Rb,B)[1],min_EbN0_mod(BER,mod)[1])
    ebn0 = EbN0(Pt,Gt,Grlin,f,hkm,elev,Ta,FdB,Rb,LldB,RE) 
    marg = ebn0[1] - minEbN0
    print('minEbN0(dB)= {}, link margin (dB) = {}'.format(minEbN0,marg))
    return marg

def size_antenna(Pt,Gr,f,hkm,elev,Ta,FdB,Rb,B,BER,margindB=0,mod = 'QPSK',LldB=0,RE=R_earth):
    from scipy.optimize import minimize_scalar
    err = lambda Dt:(10000*np.abs( link_margin(Pt,Dt,Gr,f,hkm,elev,Ta,FdB,Rb,B,BER,mod = 'QPSK',LldB=0,RE=R_earth)- margindB))
    res = minimize_scalar(err, bounds = (0,10))
    Dt_opt = res.x
    return Dt_opt