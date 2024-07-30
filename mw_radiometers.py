# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 17:33:20 2018

@author: Dani Selva
"""
from numpy import cos,tan,radians,pi,arccos,sin,abs
from orbits import R_earth


## Swath
def flat_swath(h,eta):
    
    sw = 2*h*tan(radians(eta))
    return sw

def curved_swath(h,eta,RE=R_earth):
    
    epsilon = arccos(sin(radians(eta))*(RE+h)/RE)
    lamda = pi/2 - radians(eta) - epsilon
    sw = 2*RE*lamda
    return sw

def compare_swath():
    import matplotlib.pyplot as plt
    from numpy import linspace
    h = linspace(400,800, num = 50)
    etas = linspace(10,50,num=5)
    for eta in etas:
        flat_swaths = flat_swath(h,eta)
        curved_swaths = curved_swath(h,eta)
        errs = abs(flat_swaths - curved_swaths)
        plt.plot(h,errs, label = 'eta = ' + str(round(eta)))
    plt.legend()
    plt.xlabel('h(km)')
    plt.ylabel('|flat swath - curved swath|(km)')
    return max(errs)

## Spatial resolution
    
def spat_res_nadir(h,landa,D,alfa=1):
    dtheta = alfa*landa/D
    dx = h*dtheta
    return dx

def cross_track_spat_res_off_nadir(h,f,D,eta,alfa=1):
    from comms import f2l
    landa = f2l(f)
    dtheta = alfa*landa/D
    dx = h*(tan(radians(eta) + dtheta/2) - tan(radians(eta) - dtheta/2))
    return dx


def along_track_spat_res_off_nadir(h,f,D,eta,alfa=1):
    from comms import f2l
    landa = f2l(f)
    dtheta = alfa*landa/D
    dx = h/cos(radians(eta))*dtheta
    return dx


def cmis_hsr(D,eta,f,h):
    alfa = 1.22
    return cross_track_spat_res_off_nadir(h,f,D,eta,alfa)

def size_cmis(hsr,eta,f,h):
    """Returns the aperture of a MWR required to achieve the desired hsr for a given off-nadir angle eta, altitude h, and frequency f"""
    from scipy.optimize import minimize
    
    err = lambda D:(abs(cmis_hsr(D,eta,f,h) - hsr))
    D0=1
    res = minimize(err,D0)
    return res.x    


## received power and SNR
def Teq(f,D,T,e_avg):
    from comms import f2l
    l=f2l(f)
    A=pi*D**2/4 # circular aperture
    return T*A/(l**2)*e_avg

def received_power(Teq,B):
    from remote_sensing import kb
    return kb*Teq*B
