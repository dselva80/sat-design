# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 10:32:59 2018

@author: dselva
"""
import numpy as np
from orbits import R_earth,mu_earth

def Fnum(f,D):
    return f/D

def diff_HSR(landa,h,D,eta=0,circ=True):
    """Returns cross-track,along-track diffraction limited HSR for circular aperture telescope"""
    if circ:
        alpha = 1.22
    else:
        alpha = 1.0
        
    hsr_nadir = alpha*landa*h/D
    return hsr_nadir/np.cos(np.radians(eta)),hsr_nadir/(np.cos(np.radians(eta))**2)

def GSD(p,f,h):
    return p*h/f
    
def Qfactor(f,p,h,D,landa,eta=0):
    return diff_HSR(landa,h,D,eta)[0]/GSD(p,f,h)

def data_volume_per_scan_GEO(pix_size,Nband,bit_per_pixel):
    "Returns the number of bits for one full scan of the Earth surface giventhe pixel size number of bands and bits per pixel"
    from geometry import cone_solid_angle,planet_subtended_angle
    IFOV = np.degrees(2*np.arctan(pix_size/36000000))
    Npixel = cone_solid_angle(planet_subtended_angle(36000,6370))/cone_solid_angle(IFOV)
    DV = Npixel*Nband*bit_per_pixel
    return DV
    
def data_rate_pushbroom(pix_size,Nband,bit_per_pixel,h,Npix,RE=R_earth,mu=mu_earth):
    from orbits import orb_vel_c
    V=orb_vel_c(h,mu,RE)
    Vg=V*RE/(RE+h)
    dt = pix_size/Vg
    R= Nband*Npix*bit_per_pixel/dt
    return R

def data_rate_whiskbroom(pix_size,scan_angle,Nband,bit_per_pixel,h,RE=R_earth,mu=mu_earth):
    from orbits import orb_vel_c
    IFOV= 2*np.arctan(pix_size/2/h/1000)
    Npix = np.radians(scan_angle)/IFOV
    V=orb_vel_c(h,mu,RE)
    Vg=V*RE/(RE+h)
    dt = pix_size/Vg
    R= Nband*Npix*bit_per_pixel/dt
    return R