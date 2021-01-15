#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 07-12-2020
# alex
# test_relE.py

import sympy as sp
from sympy.physics.quantum import Operator, Commutator
import sys
sys.path.insert(0, '../')
import relations as rl
from qolib import *

if __name__ == "__main__":
    A = Operator("A")
    B = Operator("B")

    E = llegir_base("niats_base.txt", A, B)
    ord_bch = len(rl.tamE) + 1
    for i in range(1, ord_bch):
        for j in range(1, rl.tamE[i - 1] + 1):
            for k in range(1, ord_bch):
                for l in range(1, rl.tamE[k - 1] + 1):
                    if rl.relE[i][j][k][l][0][0] != 0:
                        res = rl.relE[i][j][k][l]
                        cad = "[E" + str(i) + str(j) + ", E" + str(k) + str(l) + "] - ("
                        primer = True
                        resq = sp.S(0)
                        for m in res:
                            resq += sp.Rational(m[0], m[1])*E[m[2]][m[3]]
                            if m[0] > 0 and primer == False:
                                cad += "+ "
                            elif m[0] < 0:
                                cad += "- "
                            primer = False
                            if abs(m[0]) > 1:
                                cad += str(abs(m[0])) + "*"
                            cad += "E" + str(m[2]) + str(m[3]) + ""
                            if m[1] != 1:
                                cad += "/" + str(m[1])
                            cad += " "
                        cad = cad[:-1] + ") ="
                        r = (Commutator(E[i][j], E[k][l]) - resq).doit().expand()
                        if (r != 0):
                            print(cad, r)
