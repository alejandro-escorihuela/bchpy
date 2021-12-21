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

def baserkn9(w):
    return w[1:10] + w[11:27]

if __name__ == "__main__":
    a = sp.symbols("a0:14")
    b = sp.symbols("b0:15")

    printi("Calculant mètode RKNA19")
    t0 = tm.time()
    na19 = Metode(9)
    na19.setABA(a[0], b[0], a[1], b[1], a[2], b[2], a[3], b[3], a[4], b[4], a[5], b[5], a[6], b[6], a[7], b[7], a[8], b[8], a[9], b[9],
               a[9], b[8], a[8], b[7], a[7], b[6], a[6], b[5], a[5], b[4], a[4], b[3], a[3], b[2], a[2], b[1], a[1], b[0], a[0], debug=True)  
    el9 = baserkn9(na19.w[9])
    n2 = sp.S(0)
    for i in range(len(el9)):
        printi("Calculant n2: " + str(i))
        n2 += el9[i]**2
    printi("Convertint-ho a C")
    #strC = str(sp.ccode(n2)).replace("pow", "POTENCIA")
    strC = str(n2) ## Inviable!!!
    printi("Guardant fitxer C")
    with open("efC.c", 'w') as fC:
        fC.write("#include <stdio.h>\n")
        fC.write("#include <stdlib.h>\n\n")
        fC.write("double ef9(double * a, double * br, double * bi);\n\n")
        fC.write("double ef9(double * a, double * br, double * bi) {\n")
        fC.write("  return " + strC + ";\n}\n")
        
    # t1 = tm.time()
    # print(t1 - t0, "s")
    # printi("Exportant equacions a mathematica")
    # na19.exportC("na19.c")
    # t2 = tm.time()
    # print(t2 - t1, "s")   

