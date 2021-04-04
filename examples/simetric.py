#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 09-12-2020
# alex
# simetric.py

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

    printi("Calculant mètode simètric")
    t0 = tm.time()
    met = Metode(7)
 
    #met.setABA(a[0], b[0], a[1], b[1], a[2], b[2], a[3], b[3], a[4], sp.Rational(1, 2) - (b[0] + b[1] + b[2] + b[3]), 1 - 2*(a[0] + a[1] + a[2] + a[3] + a[4]), sp.Rational(1, 2) - (b[0] + b[1] + b[2] + b[3]), a[4], b[3], a[3], b[2], a[2], b[1], a[1], b[0], a[0], debug=True)
    # a13, b13 = sp.symbols("a13 b13")
    # a13 = sp.Rational(1, 2) - (a[0] + a[1] + a[2] + a[3] + a[4] + a[5] + a[6] + a[7] + a[8] + a[9] + a[10] + a[11] + a[12])
    # b13 = 1 - 2*(b[0] + b[1] + b[2] + b[3] + b[4] + b[5] + b[6] + b[7] + b[8] + b[9] + b[10] + b[11] + b[12])
    met.setABA(a[0], b[0], a[1], b[1], a[2], b[2], a[3], b[3], a[4], b[4], a[5], b[5], a[6], b[6], a[7], b[7], a[8], b[8], a[9], b[9], a[10], b[10], a[11], b[11], a[12], b[12], a[13], b[13], a[13], b[12], a[12], b[11], a[11], b[10], a[10], b[9], a[9], b[8], a[8], b[7], a[7], b[6], a[6], b[5], a[5], b[4], a[4], b[3], a[3], b[2], a[2], b[1], a[1], b[0], a[0], debug=True)  
    t1 = tm.time()
    print(t1 - t0, "s")

    printi("Expandint mètode")
    #met.expand(debug = True)
    t2 = tm.time()
    print(t2 - t1, "s")

    # printi("Desant fitxer amb el mètode")
    # met.save("sim6_10.dat")
    # t3 = tm.time()
    # print(t3 - t2, "s")

    printi("Exportant equacions a python")
    met.exportpy("sim8_27.py")
    t2 = tm.time()
    print(t2 - t1, "s")
    

