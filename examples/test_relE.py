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
    
    for i in range(len(rl.tamE)):
        for j in range(rl.tamE[i - 1]):
            for k in range(len(rl.tamE)):
                for l in range(rl.tamE[k]):
                    if rl.relE[i + 1][j + 1][k + 1][l + 1][0][0] != 0:
                        res = rl.relE[i + 1][j + 1][k + 1][l + 1]
                        cad = "[E" + str(i + 1) + str(j + 1) + ", E" + str(k + 1) + str(l + 1) + "] - ("
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
                        r = (Commutator(E[i + 1][j + 1], E[k + 1][l + 1]) - resq).doit().expand()
                        print(cad, r)
