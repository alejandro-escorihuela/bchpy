#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 13-01-2021
# alex
# bchterms.py

import sympy as sp
from sympy.physics.quantum import Operator, Commutator
import time as tm
import sys
sys.path.insert(0, '../')
from qolib import *

if __name__ == "__main__":
    n = 10
    A = Operator("A")
    B = Operator("B")
    t0 = tm.time()
    bchexpr = bchOp(A, B, n)
    t1 = tm.time()
    print(bchexpr)
    print(t1 - t0, "s")
    bchvec = [sp.S(0)]*n
    for arg in bchexpr.args:
        bchvec[countDim(arg, A, B)] += arg

    for i in range(0, len(bchvec)):
        print(i, "->", bchvec[i])


