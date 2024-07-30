# -*- coding: utf-8 -*-
"""
Created on Thu May 24 09:37:07 2018

@author: Dani Selva
"""
from orbits import R_earth, mu_earth, h2T
import numpy as np
sigma = 5.67e-8
alpha_VNIRs = {'Al':0.09,'white-paint':0.2,'black-paint':0.92,'Si-Teflon':0.08,'Al-Kapton':0.38}
alpha_TIRs = {'Al':0.03,'white-paint':0.92,'black-paint':0.89,'Si-Teflon':0.8,'Al-Kapton':0.67}
dec_sun = 23.5
phi_s = 1368


#def beta()

def beta_star(h,R=R_earth):
    """"Beta Angle at which eclipses start"""
    from geometry import planet_subtended_angle
    return planet_subtended_angle(h,R)

def sun_vector(sun_azim, sun_elev):
    """Sun vector in the Earth-centric inertial frame"""
    return np.array([np.cos(np.radians(sun_azim)), np.sin(np.radians(sun_azim))*np.cos(np.radians(sun_elev)), np.sin(np.radians(sun_azim))*np.sin(np.radians(sun_elev))])

def orbit_vector(raan, i):
    return np.array([np.sin(np.radians(raan))*np.sin(np.radians(i)), -np.cos(np.radians(raan))*np.sin(np.radians(i)), np.cos(np.radians(i))])

def beta(sun_azim, sun_elev, raan, i):
    s = sun_vector(sun_azim, sun_elev)
    o = orbit_vector(raan, i)
    cos_phi = o.dot(s)
    phi = np.arccos(cos_phi)
    beta = np.degrees(phi - np.pi/2 )
    return beta

def frac_sunlight(beta, h, mu=mu_earth,R=R_earth):
    sin_theta = np.sqrt( ((R/(R+h))**2 - (np.sin(np.radians(beta)))**2) / (np.cos(np.radians(beta))**2) )
    theta = np.arcsin(sin_theta)
    T = h2T(h,mu,R)
    alfa = 2*(np.pi - theta)/(2*np.pi) # fraction of sunlight
    ecl_t = (1-alfa)*T
    return alfa,ecl_t
     
def stefan_boltzmann(T,A,e=1.0):
    """Computes the total radiant energy in W radiated by a surface A at temperature T with an emissivity e"""
    Q = sigma*e*T**4*A
    return Q

def qsun(AU_sun,alphaVNIR,A,S0=phi_s):
    """Computes the total absorbed radiant energy from the Sun on a surface A with absorptivity alpha in the VNIR at a distance AU_sun [AU] from the sun"""
    phi = S0/(AU_sun**2)
    q = phi*alphaVNIR*A
    return q

def q_alb(AU_sun,alb,alphaVNIR,A,S0=phi_s):
    """"Computes the total absorbed  solar radiance reflected off a surface with albedo alb onto a surface A with absorptivity alpha in the VNIR"""
    return alb*qsun(AU_sun,alphaVNIR,A,S0)
    
def qrad(eps, T, A):
    """Calculates the total outgoing radiant energy from a surface A with emissivity eps at temperature T"""
    return stefan_boltzmann(T,A,eps)

def vf_sphere_point(Rp,h):
    vf = (Rp/(Rp+h))**2
    return vf

def qtherm(alphaIR, Tp, ep, A, Rp, h):
    """"Computes the total absorbed radiant energy from thermal emission of a planet of radius Rp at temperature Tp with emissivity in the TIR ep, by a surface A with absorptivity alp in the TIR at an altitude h"""
    return alphaIR*vf_sphere_point(Rp,h)*ep*sigma*Tp**4*A


def dq(T,Asun, AEarth, Acold, alphaIR, alphaVNIR, eps_TIR, ep, Tp, Rp, h, AU_sun,alb,qint):
    qin = qsun(AU_sun,alphaVNIR,Asun) + q_alb(AU_sun,alb,alphaVNIR,AEarth) + qtherm(alphaIR, Tp, ep, AEarth, Rp, h) + qint
    qout= qrad(eps_TIR, T, Acold)
    dq = abs(qout - qin)
    return dq


def eq_temp_simple(a,eps,As,Ar,AU=1):
    T4 = phi_s/sigma * a/eps * As/Ar
    return T4**0.25
    
def eq_temp(Asun, AEarth, Acold, alphaIR, alphaVNIR, eps_TIR, ep, Tp, Rp, h, AU_sun,alb,qint):
    from scipy.optimize import minimize
    fun = lambda T:dq(T,Asun, AEarth, Acold, alphaIR, alphaVNIR, eps_TIR, ep, Tp, Rp, h, AU_sun,alb,qint)
    t0 = 300
    res = minimize(fun, t0, method='nelder-mead',options={'xtol': 1e-3, 'disp': True})
    return res.x

def rad_cool_t(m,Ti,Tf,e,A,cp=900):
    dt = m*cp*(Tf**-3-Ti**-3)/3/sigma/A/e
    return dt

def rad_cool_temp(m,Ti,dt,e,A,cp=900):
    Tf3 = Ti**-3 + e*sigma*A*dt/m/cp
    return (Tf3)**(-1.0/3)