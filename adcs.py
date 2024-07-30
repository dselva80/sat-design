# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 11:16:28 2018

@author: dselva
"""
from numpy import radians,abs,sin,cos,sqrt

G = 6.67e-11
M = 5.98e24
RE = 6378
    
def deg2arcsec(deg):
    return deg*3600

def deg2microrad(deg):
    return radians(deg)*1e6

def disturbance_torques(h,Cd,Adrag,cp_cg,D,Iz,Iy,off_nadir,q,cps_cg,Asun,sun_angle,mag_colat=0):
    """Computes all disturbance torques and the worst case"""
    tau_drag = drag_dist_torque(h,Cd,Adrag  ,cp_cg)
    tau_mg = magn_dist_torque(h,D,mag_colat)
    tau_grav = grav_dist_torque(h,Iz,Iy,off_nadir)
    tau_sun = solar_torque(q,cps_cg,Asun,sun_angle)
    torque = max(tau_drag,tau_mg,tau_grav,tau_sun)
    torques = {"grav":tau_grav,"magn":tau_mg,"aero":tau_drag,"sun":tau_sun}
    return torque,torques


def drag_dist_torque(h,Cd,A,cp_cg):
    from atmosphere import drag
    tau = drag(h,Cd,A)*cp_cg
    return tau

def magn_dist_torque(h,D,mag_colat=0):
    R = 1000*(RE + h)
    M = 7.8e15
    lamda = sqrt(1+3*cos(radians(mag_colat))**2)
    B = lamda*M/(R**3)
    tau = D*B
    return tau

def grav_dist_torque(h,Iz,Iy,off_nadir):
    r = 1000*(RE + h)
    mu = G*M
    tau = 1.5*(mu/(r**3))*abs(Iz-Iy)*sin(radians(2*off_nadir))
    return tau

def solar_torque(q,cps_cg,A,sun_angle):
    c = 3e8
    phi_s = 1370
    F = phi_s/c*A*(1+q)*cos(radians(sun_angle))
    tau = F*cps_cg
    return tau

def inertia_cuboid(a,b,c,msat,d=0):
    import numpy as np
    Isat_xx = msat/12*(b**2 + c**2)
    Isat_yy = msat/12*(a**2 + c**2)
    Isat_zz = msat/12*(a**2 + b**2)
    Isat = np.array([Isat_xx, Isat_yy, Isat_zz])
    return Isat

def inertia_sat_solar_panel(a,b,c,msat,d,e,f,t,msol):
    "Calculates the mass moment of inertia of a cuboid satellite of dimensions a,b,c with 2 symmetric solar panels at distance d from the center of mass of the cuboid and with dimensions e,f,t"
    import numpy as np
    Isat_xx = msat/12*(b**2 + c**2)
    Isat_yy = msat/12*(a**2 + c**2)
    Isat_zz = msat/12*(a**2 + b**2)
    Isat = np.array([Isat_xx, Isat_yy, Isat_zz])
    
    Isol_cm_xx = msol/12*(f**2 + t**2)
    Isol_cm_yy = msol/12*(e**2 + t**2)
    Isol_cm_zz = msol/12*(e**2 + f**2)
    Isol_cm = np.array([Isol_cm_xx,Isol_cm_yy,Isol_cm_zz])
    Isol = Isol_cm + np.array([msol*d**2, 0, msol*d**2])
    
    I = Isat  + 2*Isol
    return I

def slewing_torque(I,dtheta,dt):
   "Calculate the torque required to slew assuming a triangular profile for omega"    
   alpha = 4*radians(dtheta)/(dt**2)
   torque = I*alpha
   return torque
   
def momentum_storage(T,h):
    from orbits import h2T
    P = h2T(h)
    h = T*P/4/sqrt(2)
    return h