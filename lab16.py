# -*- coding: utf-8 -*-
"""
Created on Sun Nov 11 16:56:06 2018

@author: Ana-Dani
"""
import numpy as np
def mbool(b):
    if b == '1': 
        return True
    else: 
        return False
        
def switch_f(state):
    ## state = [s1,s2,c1,a1,a2,a3]
    s1 = mbool(state[0])
    s2 = mbool(state[1])
    c1 = mbool(state[2])
    a1 = mbool(state[3])
    a2 = mbool(state[4])
    a3 = mbool(state[5])
    X = 1- (1-s1*s2)*(1-c1)*(1-a1*a2*a3)
    return X

def de2bi(i,N):
    s = '\'{0:0' + str(N) + 'b}\'.format(i)'
    b = eval(s)
    return b

def state_prob(x,vec_r):
    N = np.size(vec_r)
    p = 1.0
    for i in range(N):
        if x[i] == '1':
            p = p*(1-vec_r[i])
        else:
            p = p*(vec_r[i])    
    return p
        
def reliab_bf(vec_r,sf):
    """Calculates the reliability of a system given a vector of failure rates and a switching function"""
    N = np.size(vec_r)
    R = 0.0
    for d in range(2**N):
        x = de2bi(d,N)
        fail = switch_f(x)
        if not fail:
            p = state_prob(x,vec_r)
            R = R + p
    
    return R