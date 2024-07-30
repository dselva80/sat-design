# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 12:01:18 2018

@author: dselva
"""
from orbits import R_earth,mu_earth,R_titan,mu_titan,orb_vel_c,hkm2a
import numpy as np

def isothermal_density(h,rho_ref=4.13e-13,h_ref=500,H=41.67):
    from numpy import exp
    rho = rho_ref*exp(-(h-h_ref)/H)
    return rho

def drag(h,Cd,A,rho_ref=4.13e-13,h_ref=500,H=41.67):
    V = orb_vel_c(h)
    rho = isothermal_density(h,rho_ref=4.13e-13,h_ref=500,H=41.67)
    D = 0.5*rho*V*V*A*Cd
    return D

def alt_loss_1rev(h0,Cd,ballistic,R=R_earth,rho_ref=4.13e-13,h_ref=500,H=41.67):
    a = hkm2a(h0,R)
    rho = isothermal_density(h0,rho_ref,h_ref,H)
    da = 2*np.pi*Cd/ballistic*a**2*rho
    return da
    