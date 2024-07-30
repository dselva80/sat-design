# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 22:12:01 2018

@author: Dani Selva
"""

def learning_curve(S,N,C1):
    from numpy import log
    B = 1+log(S)/log(2)
    Cn = C1*N**B
    
    return Cn,Cn/N

def test_learn_curv(S,N,C1):
    from numpy import log
    B = 1+log(S)/log(2)
    cs = [C1]
    cumcost = C1
    for i in range(2,N+1):
        Ci = C1*i**B
        ci = Ci-cumcost
        cs.append(ci)
        cumcost+= ci
    return cs
        
    