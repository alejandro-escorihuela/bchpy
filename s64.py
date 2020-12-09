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
    print("BCH")
    x, h = sp.symbols("x h")
    a = sp.symbols("a0:3")
    b = sp.symbols("b0:4")
    t0 = tm.time()
    met = schemeBAB(b[0], a[0], b[1], a[1], b[2], sp.Rational(1, 2) - a[1] - a[0], sp.S(1) - (sp.S(2)*(b[0] + b[1] + b[2])), sp.Rational(1, 2) - a[1] - a[0], b[2], a[1], b[1], a[0], b[0])
    w = mat_sch(met)
    for i in range(1, len(w)):
        for j in range(1, len(w[i])):
            print("w(" + str(i) + "," + str(j) + ") =", w[i][j])
    t1 = tm.time()
    print(t1-t0, "s")
