# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 10:08:29 2018

@author: dselva
"""
import numpy as np
from orbits import R_earth, mu_earth, hkm2a, J2_earth
Omegadot_sun = 2*np.pi/365.25/86400


def Omegadot(a,i,e,mu=mu_earth,R=R_earth,J2=J2_earth):# RAAN dot
    from orbits import a2n
    n = a2n(a,mu)
    odot = -1.5*n*J2*((1000*R/a)**2)*np.cos(np.radians(i))/((1-e**2)**2)
    return odot

def SSOh2i(hkm,e=0,mu=mu_earth, R=R_earth, J2=J2_earth):
    from scipy.optimize import minimize
    err = lambda i:(1e10*abs(Omegadot(hkm2a(hkm,R),i,e,mu,R,J2) - Omegadot_sun))
    i0=90
    res = minimize(err,i0)
    return res.x[0]

def omegadot(a,i,e,mu=mu_earth,R=R_earth,J2=J2_earth):# arg perigee dot
    from orbits import a2n
    n = a2n(a,mu)
    odot = 0.75*n*J2*((1000*R/a)**2)*(5*np.cos(np.radians(i))**2-1)/((1-e**2)**2)
    return odot

def deltan(a,i,e,mu=mu_earth,R=R_earth,J2=J2_earth):# ndot - n0
    from orbits import a2n
    n = a2n(a,mu)
    dn = 0.75*n*J2*((1000*R/a)**2)*(3*np.cos(np.radians(i))**2-1)/((1-e**2)**1.5)
    return dn
