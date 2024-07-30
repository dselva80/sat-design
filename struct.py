# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 23:32:45 2018

@author: Ana-Dani
"""
import numpy as np

def bend_mode(E,I,m,L):
    k = 3*E*I/L**3
    omega = np.sqrt(k/m)
    f = omega/2/np.pi
    return f

def area_moment_ring(do,di):
    return np.pi/64*(do**4 - di**4)

def tors_stiff_boom(G,r,t,L):
    J = 2.0/3*np.pi*r*t**3
    return G*J/L
    
def solar_panel_mode(k,m,a,b,d):
    from numpy import sqrt,pi
    J0=m*(a**2+b**2)/12    
    J = J0+m*d**2    
    omega = sqrt(k/J)
    f = omega/2/pi
    return f

def buckling_load(E,I,k,L):
    Lp = k*L
    P = np.pi**2*E*I/Lp**2
    return P

def cad():
    import pandas as pd
    from shutil import copyfile

    
    df = pd.read_excel("cad.xlsx")
    
    #n = len(df.index) - 1
    copyfile("cylshell.txt", "cad.txt")
    fi = open("cad.txt","a+")
    M = 0
    xcm = 0
    ycm = 0
    zcm = 0
    for index, row in df.iterrows():

        x = row['xc']
        y = row['yc']
        z = row['zc']
        ty = row['type']
        m = row['m']
        M = M+m
        xcm = xcm + x*m
        ycm = ycm + y*m
        zcm = zcm + z*m
        
        if ty == "cube":
            a = row['a']
            b = row['b']
            c = row['c']
            str = "translate([{},{},{}]) cube([{},{},{}]);\n".format(x,y,z,a,b,c)
        elif ty == "cylinder":
            r = row['r']
            h = row['h']
            str = "translate([{},{},{}]) cylinder(r={},h={},$fn=20);\n".format(x,y,z,r,h)
        elif ty == "sphere":
            r = row['r']
            str = "translate([{},{},{}]) sphere({},$fn=20);\n".format(x,y,z,r)
        elif ty == "cylshell":
            r = row['r']
            h = row['h']
            t = row['t']
            str = "translate([{},{},{}]) cylshell({},{},{});\n".format(x,y,z,r,h,t)
        elif ty == "sphershell":
            r = row['r']
            t = row['t']
            str = "translate([{},{},{}]) sphershell({},{});\n".format(x,y,z,r,t)
        fi.write(str)
    xcm = xcm / M
    ycm = ycm / M
    zcm = zcm / M
    fi.close()
    return xcm,ycm,zcm