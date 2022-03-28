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
from bchpy import printi

if __name__ == "__main__":
    n = 8
    A = Operator("A")
    B = Operator("B")

    printi("Expansió BCH")
    t0 = tm.time()
    bchexpr = bchOp(A, B, n)
    t1 = tm.time()
    print(bchexpr)
    print(t1 - t0, "s")
    bchvec = [sp.S(0)]*n
    for arg in bchexpr.args:
        bchvec[countDim(arg, A, B)] += arg        
    E = llegir_base("niats_base.txt", A, B)
    print()
    
    printi("BCH en termes dels elements de la base")
    t0 = tm.time()
    for i in range(1, len(bchvec)):
        c = sp.symbols("s1:" + str(len(E[i])))
        cesq = sp.S(0)
        for m in range(1, len(E[i])):
            cesq += c[m - 1]*E[i][m]
        cesq = cesq.doit().expand()
        sol = resoldre(cesq, bchvec[i])
        cad = ""
        nt = 0
        for m in sol:
            if sol[m] != 0:
                nume, deno = sp.fraction(sol[m])
                frac = sp.Rational(nume, deno)
                if frac > 0:
                    cad += "+"
                cad += str(frac) + "*E[" + str(i) + "][" + str(c.index(m) + 1) + "]"
                # cad += str(frac) + "*" + str(E[i][c.index(m) + 1])
                nt += 1
        print(str(i) + "(" + str(nt) + ") ->", cad)
    t1 = tm.time()
    print(t1 - t0, "s")
