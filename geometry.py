# -*- coding: utf-8 -*-
"""
Created on Sat Mar 24 20:12:59 2018

@author: Dani Selva
"""
import numpy as np
from numpy import sin,cos,arctan,arcsin,sqrt,arccos,degrees,radians

def planet_subtended_angle(h,RE):
    """Calculates the angle subtended by the Earth from the spacecraft"""
    sinrho =RE/(RE+h)
    rho = degrees(arcsin(sinrho))
    return rho

def azel2xyz(RE,azdeg,eldeg):
    az = radians(azdeg)
    el = radians(eldeg)
    x = RE*cos(az)*cos(el)
    y = RE*sin(az)*cos(el)
    z = RE*sin(el)
    return x,y,z

def xyz2azel(x,y,z):
    RE = sqrt(x**2 + y**2 + z**2)
    if x>0:
        az = degrees(arctan(y/x))
    else:
        az = 90    
    el = degrees(arcsin(z/RE))
    if el==90.0:
        az = np.nan
        
    return az,el

def eclipse_time(h,RE,beta):
    """Calculates the duration of eclipse for a given orbit and sun angle condition"""
    from orbits import h2T
    T = h2T(h)
    rho = planet_subtended_angle(h,RE)
    cos_phi2 = cos(radians(rho))/cos(radians(beta))
    phi = 2*arccos(cos_phi2)
    f = phi/(2*np.pi)
    ecl_t = f*T
    return ecl_t,f
    
def project_fov(half_fov,h,RE):
    "Calculates the Earth central angle lambda [deg] and ground swath deltax in the same units as h, using the law of sines"
    sin_rho =RE/(RE + h)
    cos_eps = sin(radians(half_fov))/sin_rho
    eps= degrees(arccos(cos_eps))
    #eta = arcsin(cos_eps*sin_rho)
    lamda = 90 - half_fov - eps
    deltax = 2*RE*radians(lamda)
    return lamda,deltax

def incidence_angle(off_nadir,h,RE):
    sin_rho =RE/(RE + h)
    cos_eps = sin(radians(off_nadir))/sin_rho
    eps= degrees(arccos(cos_eps))
    return eps

def eta(eps,RE,h):
    "Calculates the look angle or off nadir angle eta (deg) from a spacecraft that corresponds to a given altitude and elevation or incidence angle from the ground eps [deg]"
    sin_eta = RE/(RE+h)*np.cos(np.radians(eps))
    return np.degrees(np.arcsin(sin_eta))
    
def comm_range(hkm,eta,RE):
    """ Calculates range given min elev angle eta, orbit altitude h, and planet radius"""
    if eta == 90:
        return hkm
    sin_eps = RE/(RE+hkm)*np.cos(np.radians(eta))
    eps = np.arcsin(sin_eps)
    sin_alpha= np.cos(radians(eta)+eps)
    R = RE*sin_alpha/sin_eps
    return R,degrees(eps),degrees(np.arcsin(sin_alpha))

def cone_solid_angle(theta):
    return 2*np.pi*(1-np.cos(np.radians(theta)))