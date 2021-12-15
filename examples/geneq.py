#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 05-02-2021
# alex
# geneq17.py

import numpy as np
import sympy as sp
import time as tm
import sys
sys.path.insert(0, '../')
from bchpy import *
import relations as rl

if __name__ == "__main__":
    a = sp.symbols("a0:14")
    b = sp.symbols("b0:15")
  
    printi("Calculant mètode RKNA17")
    t0 = tm.time()
    na17 = Metode(7)
    na17.setABA(a[0], b[0], a[1], b[1], a[2], b[2], a[3], b[3], a[4], b[4], a[5], b[5], a[6], b[6], a[7], b[7], a[8], b[8],
                a[8], b[7], a[7], b[6], a[6], b[5], a[5], b[4], a[4], b[3], a[3], b[2], a[2], b[1], a[1], b[0], a[0], debug=True)  
    t1 = tm.time()
    print(t1 - t0, "s")
    printi("Exportant equacions a mathematica")
    na17.exportnb("na17")
    t2 = tm.time()
    print(t2 - t1, "s")
    
    printi("Calculant mètode RKNA18")
    t0 = tm.time()
    na18 = Metode(7)
    na18.setABA(a[0], b[0], a[1], b[1], a[2], b[2], a[3], b[3], a[4], b[4], a[5], b[5], a[6], b[6], a[7], b[7], a[8], b[8], a[9],
                b[8], a[8], b[7], a[7], b[6], a[6], b[5], a[5], b[4], a[4], b[3], a[3], b[2], a[2], b[1], a[1], b[0], a[0], debug=True)  
    t1 = tm.time()
    print(t1 - t0, "s")
    printi("Exportant equacions a mathematica")
    na18.exportnb("na18")
    t2 = tm.time()
    print(t2 - t1, "s")

    printi("Calculant mètode RKNA19")
    t0 = tm.time()
    na19 = Metode(7)
    na19.setABA(a[0], b[0], a[1], b[1], a[2], b[2], a[3], b[3], a[4], b[4], a[5], b[5], a[6], b[6], a[7], b[7], a[8], b[8], a[9], b[9],
                a[9], b[8], a[8], b[7], a[7], b[6], a[6], b[5], a[5], b[4], a[4], b[3], a[3], b[2], a[2], b[1], a[1], b[0], a[0], debug=True)  
    t1 = tm.time()
    print(t1 - t0, "s")
    printi("Exportant equacions a mathematica")
    na19.exportnb("na19")
    t2 = tm.time()
    print(t2 - t1, "s")    

    printi("Calculant mètode RKNB17")
    t0 = tm.time()
    nb17 = Metode(7)
    nb17.setBAB(b[0], a[0], b[1], a[1], b[2], a[2], b[3], a[3], b[4], a[4], b[5], a[5], b[6], a[6], b[7], a[7], b[8], a[8],
                b[8], a[7], b[7], a[6], b[6], a[5], b[5], a[4], b[4], a[3], b[3], a[2], b[2], a[1], b[1], a[0], b[0], debug=True)  
    t1 = tm.time()
    print(t1 - t0, "s")
    printi("Exportant equacions a mathematica")
    nb17.exportnb("nb17")
    t2 = tm.time()
    print(t2 - t1, "s")    

    printi("Calculant mètode RKNB18")
    t0 = tm.time()
    nb18 = Metode(7)
    nb18.setBAB(b[0], a[0], b[1], a[1], b[2], a[2], b[3], a[3], b[4], a[4], b[5], a[5], b[6], a[6], b[7], a[7], b[8], a[8], b[9],
                a[8], b[8], a[7], b[7], a[6], b[6], a[5], b[5], a[4], b[4], a[3], b[3], a[2], b[2], a[1], b[1], a[0], b[0], debug=True)  
    t1 = tm.time()
    print(t1 - t0, "s")
    printi("Exportant equacions a mathematica")
    nb18.exportnb("nb18")
    t2 = tm.time()
    print(t2 - t1, "s")  
    
    printi("Calculant mètode RKNB19")
    t0 = tm.time()
    nb19 = Metode(7)
    nb19.setBAB(b[0], a[0], b[1], a[1], b[2], a[2], b[3], a[3], b[4], a[4], b[5], a[5], b[6], a[6], b[7], a[7], b[8], a[8], b[9], a[9],
                b[9], a[8], b[8], a[7], b[7], a[6], b[6], a[5], b[5], a[4], b[4], a[3], b[3], a[2], b[2], a[1], b[1], a[0], b[0], debug=True)  
    t1 = tm.time()
    print(t1 - t0, "s")
    printi("Exportant equacions a mathematica")
    nb19.exportnb("nb19")
    t2 = tm.time()
    print(t2 - t1, "s")    
