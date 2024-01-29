#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 10-01-2021
# alex
# recursive.py

import ctypes as ct
import numpy as np
import sympy as sp
import os
import sys
sys.path.insert(0, os.path.dirname(__file__) + "/recurpy")
import recurABpy as rABpy
import recurXXpy as rXXpy
import recurS2py as rS2py
import recurS4py as rS4py
import recurS6py as rS6py

rABc = ct.CDLL(os.path.dirname(__file__) + "/recurc/recurABc.so")
rABc.recAsim_c.argtypes = (ct.POINTER(ct.c_double), ct.c_double, ct.c_int, ct.c_int)    
rABc.recAsim_c.restype = ct.c_void_p
rABc.recBsim_c.argtypes = (ct.POINTER(ct.c_double), ct.c_double, ct.c_int, ct.c_int)    
rABc.recBsim_c.restype = ct.c_void_p

def recA(alp, bet, x, order):
    return rABpy.recA_py(alp, bet, x, order)

def recB(alp, bet, x, order):
    return rABpy.recB_py(alp, bet, x, order)

def recAsim(alp, bet, x, order, rkn):
    if not isinstance(x, complex) and not isinstance(x, sp.Basic):
        global rABc
        alpc = (ct.c_double*128)()
        k = 0
        for i in range(1, len(alp), 2):
            for j in range(1, len(alp[i])):
                alpc[k] = alp[i][j]
                k += 1
        rABc.recAsim_c(alpc, x, order, int(rkn))
        k = 0
        for i in range(1, len(alp), 2):
            for j in range(1, len(alp[i])):
                bet[i][j] = alpc[k]
                k += 1 
        return bet
    else:
        return rABpy.recAsim_py(alp, bet, x, order, rkn)

def recBsim(alp, bet, x, order, rkn):
    if not isinstance(x, complex) and not isinstance(x, sp.Basic):
        global rABc       
        alpc = (ct.c_double*128)()
        k = 0
        for i in range(1, len(alp), 2):
            for j in range(1, len(alp[i])):
                alpc[k] = alp[i][j]
                k += 1
        rABc.recBsim_c(alpc, x, order, int(rkn))
        k = 0
        for i in range(1, len(alp), 2):
            for j in range(1, len(alp[i])):
                bet[i][j] = alpc[k]
                k += 1 
        return bet
    else:
        return rABpy.recBsim_py(alp, bet, x, order, rkn)

def recXa(alp, bet, x, order):
    return rXXpy.recXa_py(alp, bet, x, order)

def recXb(alp, bet, x, order):
    return rXXpy.recXb_py(alp, bet, x, order)

def recS2(alp, bet, x, order):
    return rS2py.recS2_py(alp, bet, x, order)

def recS2sim(alp, bet, x, order):
    return rS2py.recS2sim_py(alp, bet, x, order)

def recS4(alp, bet, x, order):
    return rS4py.recS4_py(alp, bet, x, order)

def recS4sim(alp, bet, x, order):
    return rS4py.recS4sim_py(alp, bet, x, order)

def recS6sim(alp, bet, x, order):
    return rS6py.recS6sim_py(alp, bet, x, order)
