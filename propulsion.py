# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 08:28:29 2018

@author: dselva
"""

import numpy as np
import matplotlib.pyplot as plt

g0 = 9.81

def rocket_eq_plot(dV,Isp):
    #Isp=400
    Ve=g0*Isp
    #dV = 1e4
    mr = np.exp(dV/Ve) - 1
    plt.scatter(dV,mr)
    plt.xlabel('\Delta V')
    plt.ylabel('m_prop / m_dry')
    return

def rocket_eq_plot2(dV,Isp,eps):
    #Isp=400
    Ve=g0*Isp
    #dV = 1e4
    R = np.exp(dV/Ve) 
    lamda = R*eps-1/(1-R)
    plt.scatter(dV,lamda)
    plt.xlabel('\Delta V')
    plt.ylabel('\lambda')
    return

def rocket_exp(dV,Isp):
    Ve=g0*Isp
    #dV = 1e4
    mr = np.exp(dV/Ve) - 1
    return mr

def prop_mass(dV,Isp,mdry):
    #Isp=400
    mr = rocket_exp(dV,Isp)
    mprop = mr*mdry
    return mprop


def iterative_rocket_eq(dV,nmanuevers,Isp,wet_mass):
    mi = wet_mass
    mprop = 0
    for i in range(0,nmanuevers):
        mf=mi*rocket_exp(dV,Isp)
        mprop = mprop + mi - mf
        mi = mf
    return mprop
        
def propmass2dV(mi,mf,Isp):
    #Isp=400
    Ve=g0*Isp
    return Ve*np.log(mi/mf)


def thrust_req(dV,dt,Isp,mdry):
    mp = prop_mass(dV,Isp,mdry)
    print(mp)
    mdot = mp/dt
    print(mdot)
    T = mdot*g0*Isp
    return T

def burn_time(dV,T,Isp,mdry):
    mp = prop_mass(dV,Isp,mdry)
    print(mp)
    mdot = T/g0/Isp
    print(mdot)
    dt = mp/mdot
    return dt