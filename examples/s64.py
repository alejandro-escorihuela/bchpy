#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 09-12-2020
# alex
# s64.py

import numpy as np
import sympy as sp
import relations as rl
import time as tm
from bchpy import *

if __name__ == "__main__":
    print("s64")
    a = sp.symbols("a0:3")
    b = sp.symbols("b0:4")
    
    t0 = tm.time()
    met = Metode(5)
    met.setABA(b[0], a[0], b[1], a[1], b[2], sp.Rational(1, 2) - a[1] - a[0], sp.S(1) - (sp.S(2)*(b[0] + b[1] + b[2])), sp.Rational(1, 2) - a[1] - a[0], b[2], a[1], b[1], a[0], b[0], debug=True)
    # met.setABA(sp.S(0.0792036964311956), sp.S(0.2095151066133620), sp.S(0.3531729060497740), sp.S(-0.143851773179818), sp.S(-0.042065080357719), sp.S(0.434336666566456), sp.S(0.2193769557534987), sp.S(0.434336666566456), sp.S(-0.042065080357719), sp.S(-0.143851773179818), sp.S(0.3531729060497740), sp.S(0.2095151066133620), sp.S(0.0792036964311956), debug=True)
    met.cprint()
    t1 = tm.time()
    print(t1-t0, "s")    

