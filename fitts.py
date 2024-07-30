# -*- coding: utf-8 -*-
"""
Created on Mon May 30 09:51:47 2022

@author: danie
"""
import numpy as np
def fitts_law(W,D,t,l):
    exponent = -l*t*np.exp(-W/2/D)
    prob_fail = 1 - np.exp(exponent)
    return prob_fail
