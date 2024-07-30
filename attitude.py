# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 16:17:34 2021

@author: danie
"""
import numpy as np

def euler_angles2rot_matrix(phi,theta,psi,sequence):    
    if sequence == "313":
        cph = np.cos(np.radians(phi))
        sph = np.sin(np.radians(phi))
        cth = np.cos(np.radians(theta))
        sth = np.sin(np.radians(theta))
        cps = np.cos(np.radians(psi))
        sps = np.sin(np.radians(psi))
        
        A = np.array([[cps*cph-sps*cth*sph, cps*sph+sps*cth*cph, sps*sth],[-sps*cph-cps*cth*sph,-sps*sph+cps*cth*cph,cph*sth],[sth*sph,-sth*cph,cth]])
        return A
    else:
        return -1
    
        
    
def rot_matrix2euler_angles(A,sequence):
    if sequence == "313":
        theta = np.arccos(A[2,2])
        if np.sin(theta)>0:
            s = 1
        else:
            s = -1
                
        if np.sin(theta) != 0:
            phi = np.arctan2(s*A[2,0],-s*A[2,1])
        else:
            phi = 0
        if A[2,2] >= 0:
            phi_plus_psi = np.arctan2(A[0,1]-A[1,0],A[0,0]+A[1,1])
            psi = phi_plus_psi - phi 
        else:
            phi_minus_psi = np.arctan2(A[0,1]+A[1,0],A[0,0]-A[1,1])
            psi = phi - phi_minus_psi 
        return np.array([phi,theta,psi])
        

def rot_matrix2axis_angle(A):
    theta = np.arccos(0.5*(np.trace(A)-1))
    e = 0.5/np.sin(theta)*np.array([A[1,2]-A[2,1],A[2,0]-A[0,2],A[0,1]-A[1,0]])    
    
    return e,theta

def x_cross_mat(x):
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    return np.array([[0,-x3,x2],[x3,0,-x1],[-x2,x1,0]])

def quat_prod_mat(q):
    q1 = q[0]
    q2 = q[1]
    q3 = q[2]
    q4 = q[3]
    return np.array([[q4,q3,-q2,q1],[-q3,q4,q1,q2],[q2,-q1,q4,q3],[-q1,-q2,-q3,q4]])

def quat_mult(q1,q2):
    return quat_prod_mat(q1)*np.transpose(q2)

def vec_mult_quat(x,q):
    return quat_mult(np.array([x[0],x[1],x[2],0]),q)

def q_conj(q):
    return np.array([-q[0],-q[1],-q[2],q[3]])

def quat_rot_vec(x,q):
    tmp = vec_mult_quat(x,q_conj(q))
    tmp2 = quat_mult(q,tmp)
    return tmp2[0:3]