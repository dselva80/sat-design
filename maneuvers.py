# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 10:02:15 2018

@author: dselva
"""
import numpy as np
from orbits import hkm2a, orb_vel, orb_vel_c, hyperb_vel,a2T, mu_sun, mu_earth, R_earth, h2n, n2h

### Maneuvers    
def hohmann(h1,h2,mu=mu_earth,R=R_earth):
    aT = (hkm2a(h1,R)+hkm2a(h2,R))/2
    # dV1 from V(h1) to V(a=a1+a2/2,r=h1)
    dV1 = orb_vel(aT,hkm2a(h1,R),mu) - orb_vel_c(h1,mu,R)
    print("Burn 1: V1 = {} V2 = {} dV = {}".format(orb_vel_c(h1,mu,R),orb_vel(aT,hkm2a(h1,R),mu),dV1))
     
    # dV2 from V(a=a1+a2/2,r=h2) to V(h2)
    dV2 = orb_vel_c(h2,mu,R) - orb_vel(aT,hkm2a(h2,R),mu)
    print("Burn 2: V1 = {} V2 = {} dV = {}".format(orb_vel(aT,hkm2a(h2,R),mu),orb_vel_c(h2,mu,R),dV1))
    
    # dt and sum
    dt = 0.5*a2T(aT,mu)
    return dV1, dV2, abs(dV1) + abs(dV2),dt

def angle_change(h_km, di, mu=mu_earth,R=R_earth):
    Vi=orb_vel_c(h_km,mu,R)
    return 2*Vi*np.sin(np.radians(di/2))

def combined_maneuver(Vi,Vf,theta):
    dV = np.sqrt(Vi**2 + Vf**2 - 2*Vi*Vf*np.cos(np.radians(theta)))
    return dV

def interplanetary_hohmann(h1,h2,a1,a2,R1,R2,mu1,mu2):
    """"Perform an interplanetary Hohmann transfer between a circular orbit of altitude h1 
    around Planet 1 (semimajor axis a1, radius R1, grav constant mu1) to a circular orbit of altitude h2 around Planet 2 
    (semimajor axis a2, radius R2, grav constant mu2)"""
    
    aT= (a1 + a2)/2
    ## Burn 1: From V at circular orbit around departure planet to To V at pericenter of hyperbola around departure planet with 𝑉_∞= heliocentric velocity at departure position of transfer orbit - heliocentric velocity of departure planet 

    Vinf=orb_vel(aT,a1,mu_sun)-orb_vel(a1,a1,mu_sun)
    dV1 = hyperb_vel(hkm2a(h1,R1),Vinf) - orb_vel_c(h1,mu1)
    print("Burn 1: V1 = {} V2 = {} dV = {}".format(hyperb_vel(hkm2a(h1,R1),Vinf),orb_vel_c(h1,mu1),dV1))
    # Burn 2: From V at pericenter of hyperbola around target planet with 𝑉_∞= heliocentric velocity of target planet – heliocentric velocity at target position of transfer orbit To V at circular orbit around target planet
    Vinf2 = orb_vel(a2,a2,mu_sun) - orb_vel(aT,a2,mu_sun)
    dV2 = orb_vel(hkm2a(h2,R2),hkm2a(h2,R2),mu2) - hyperb_vel(hkm2a(h2,R2),Vinf2,mu2)
    print("Burn 1: Vf = {} Vi = {} dV = {}".format(orb_vel(hkm2a(h2,R2),hkm2a(h2,R2),mu2),hyperb_vel(hkm2a(h2,R2),Vinf2,mu2),dV2))
        
    dt = 0.5*a2T(aT,mu_sun)
    return dV1, dV2, abs(dV1)+abs(dV2), dt

def spread_sats(h_nom,h_phas,delta_max):
    dV1, dV2, dV_transfer, dt_transfer = hohmann(h_nom,h_phas)
    
    n1 = h2n(h_nom)
    n2 = h2n(h_phas)
    dn = n2 - n1
    dt_spread = np.radians(delta_max)/dn
    deltaV = 2*dV_transfer
    time = 2*dt_transfer + dt_spread
    
    print("n1 = {}, n2 = {}, dn = {}, dt_spread = {}, total dV = {}, total time = {}".format(n1,n2,dn,dt_spread,deltaV,time))
    return deltaV,time