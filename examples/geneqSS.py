#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 24-04-2023
# alex
# geneqSS.py

import numpy as np
import sympy as sp
import time as tm
import sys
sys.path.insert(0, '../')
from bchpy import *
import relations as rl

if __name__ == "__main__":
    a = sp.symbols("a0:14")
   
    printi("Calculant mètodes SS(S2)")
    printi("7 etapes simbòlic")
    t0 = tm.time()
    met0 = Metode(depth = 7, basis_type = "M", numeric = False)
    met0.setSS2(a[0], a[1], a[2], a[3], a[2], a[1], a[0], debug=True)  
    t1 = tm.time()
    print(t1 - t0, "s")
    met0.cprint()

    printi("7 etapes numèric")
    cofs = [0.78451361047755726381949763, 0.23557321335935813368479318, -1.17767998417887100694641568, 1.31518632068391121888424973]
    t0 = tm.time()
    met1 = Metode(depth = 7, basis_type = "M", numeric = True)
    met1.setSS2(cofs[0], cofs[1], cofs[2], cofs[3], cofs[2], cofs[1], cofs[0], debug = True)  
    t1 = tm.time()
    print(t1 - t0, "s")
    met1.cprint()    
    
    printi("15 etapes numèric")
    cofs = [0.74167036435061295344822780, -0.40910082580003159399730010, 0.19075471029623837995387626, -0.57386247111608226665638773,
            0.29906418130365592384446354, 0.33462491824529818378495798, 0.31529309239676659663205666, -0.79688793935291635291635401978884]
    cofs = cofs + cofs[-2::-1]
    t0 = tm.time()
    met2 = Metode(depth = 9, basis_type = "M", numeric = True)
    met2.setSS2(*cofs, debug=True)  
    t1 = tm.time()
    print(t1 - t0, "s")
    met2.cprint()  

    printi("Calculant mètodes SS(S4)")
    printi("3s etapes numèric")
    alpha = 1/(2 - 2**(1/5))
    beta = 1 -2*alpha
    cofs = [alpha, beta, alpha]
    t0 = tm.time()
    met3 = Metode(depth = 9, basis_type = "M4", numeric = True)
    met3.setSS4(*cofs, debug = True)  
    t1 = tm.time()
    print(t1 - t0, "s")
    met3.cprint()      
