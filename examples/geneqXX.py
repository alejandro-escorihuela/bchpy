#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 28-04-2021
# alex
# XX.py

import numpy as np
import sympy as sp
import time as tm
import sys
sys.path.insert(0, '../')
from bchpy import *
import relations as rl

if __name__ == "__main__":
    # sx64 = [0.0792036964311960, 0.1303114101821661, 0.2228614958676080, -0.3667132690474261, 0.3246481886897060, 0.10968847787675, 0.10968847787675, 0.3246481886897060, -0.3667132690474261, 0.2228614958676080, 0.1303114101821661, 0.0792036964311960]
    # a64 = [0.0792036964311957, 0.353172906049774, -0.0420650803577195, 0.21937695575349958]
    # b64 = [0.209515106613362, -0.143851773179818, 0.434336666566456]
    
    a = sp.symbols("a0:14")
    b = sp.symbols("b0:15")
    c = sp.symbols("c0:14")
    
    printi("Calculant composició de mètode i adjunt")
    
    metC = Metode(5, "Z")
    # metE = Metode(5, "E")

    # t0 = tm.time()
    metC.setXX(a[0], a[1], a[2], a[3], a[4], debug=True)
    # metC.setXX(*sx64, debug=True)
    # t1 = tm.time()

    # t2 = tm.time()
    # metE.setABA(a64[0], b64[0], a64[1], b64[1], a64[2], b64[2], a64[3], b64[2], a64[2], b64[1], a64[1], b64[0], a64[0], debug=True)
    # t3 = tm.time()

    # print("Composició:", t1 - t0, "s")
    # metC.cprint()
    # print("Escissió:", t3 - t2, "s")
    # metE.cprint()
    
    printi("Desant fitxer amb el mètode")
    metC.save("cad5.dat")

    printi("Exportant equacions a python")
    metC.exportpy("cad5.py")
    

