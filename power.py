# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 16:14:01 2018

@author: dselva
"""
import numpy as np
from orbits import h2T, R_earth
S0 = 1368 #W/m^2 

def size_SA(Pd,Pe,Td,Te,hkm,AU,beta,eta,D,t):
    Xd = 0.85
    Xe = 0.65
    Psa = ((Pd*Td/Xd)+(Pe*Te/Xe))/Td
    Id = 0.85*0.85*1 # temperature, design and assembly, shadowing
    Pbol = (S0/AU**2)*eta*Id*np.cos(np.radians(beta))
    Peol = Pbol*(1-D)**t
    Asa = Psa/Peol
    msa = Asa*2.5
    return msa,Asa,Pbol,Psa,Peol
    
def size_batt(Pe,Te,DOD,n=0.9):
    Cb = Pe*Te/DOD/n/3600 # C in Wh
    Msp = 100 #Wh kg
    mb = Cb/Msp
    return mb, Cb

def size_EPS(Pd,Pe,hkm,AU,beta,eta,D,t,DOD,n=0.9):
    period = h2T(hkm)
    alfa = 2.0/3
    Td = alfa*period
    Te = (1-alfa)*period
    
    msa,Asa,Pbol,Psa,Peol = size_SA(Pd,Pe,Td,Te,hkm,AU,beta,eta,D,t)
    mb, Cb = size_batt(Pe,Te,DOD,n=0.9)
    alpha = 0.25
    m = (mb + msa)* (1+alpha) 
    return m,msa,mb,Asa,Cb,Pbol,Psa,Peol

def beta_angle(hkm,R=R_earth):
    beta = np.arcsin(R/(R+hkm))
    return beta
    
def fraction_sunlight(hkm,R=R_earth):
    beta = beta_angle(hkm,R)
    cos_a = np.sqrt(hkm**2 + 2*hkm*R)/((R+hkm)*np.cos(beta))
    fe = np.arccos(cos_a)/np.pi
    return fe