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
    x, h = sp.symbols("x h")
    a = sp.symbols("a0:3")
    b = sp.symbols("b0:4")
    
    # t0 = tm.time()
    # met = esqBAB(b[0], a[0], b[1], a[1], b[2], sp.Rational(1, 2) - a[1] - a[0], sp.S(1) - (sp.S(2)*(b[0] + b[1] + b[2])), sp.Rational(1, 2) - a[1] - a[0], b[2], a[1], b[1], a[0], b[0], debug=True)
    # w = esq2mat(met)
    # for i in range(1, len(w)):
    #     for j in range(1, len(w[i])):
    #         print("w(" + str(i) + "," + str(j) + ") =", w[i][j])
    # t1 = tm.time()
    # print(t1-t0, "s")

    t0 = tm.time()
    w = mat_esqABA(b[0], a[0], b[1], a[1], b[2], sp.Rational(1, 2) - a[1] - a[0], sp.S(1) - (sp.S(2)*(b[0] + b[1] + b[2])), sp.Rational(1, 2) - a[1] - a[0], b[2], a[1], b[1], a[0], b[0], debug=True)
    # w = mat_esqABA(sp.S(0.0792036964311956), sp.S(0.2095151066133620), sp.S(0.3531729060497740), sp.S(-0.143851773179818), sp.S(-0.042065080357719), sp.S(0.434336666566456), sp.S(0.2193769557534987), sp.S(0.434336666566456), sp.S(-0.042065080357719), sp.S(-0.143851773179818), sp.S(0.3531729060497740), sp.S(0.2095151066133620), sp.S(0.0792036964311956), debug=True)
    for i in range(1, len(w)):
        for j in range(1, len(w[i])):
            print("w(" + str(i) + "," + str(j) + ") =", w[i][j])
    t1 = tm.time()
    print(t1-t0, "s")    
