#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 19-04-2023
# alex
# test_relM.py

import sympy as sp
from sympy.physics.quantum import Operator, Commutator
import sys
sys.path.insert(0, '../')
import relations as rl
from qolib import *

if __name__ == "__main__":
    n_max = 13
    S = []
    S.append([0])
    for i in range(n_max + 1):
        S.append(Operator("S_{%d}" % (i + 1)))
        
    Z = llegir_baseS2("../aux/compS2_base.txt", S)
    
    ord_bch = len(rl.tamM) + 1
    for i in range(1, ord_bch):
        for j in range(1, rl.tamM[i - 1] + 1):
            for k in range(1, ord_bch):
                for l in range(1, rl.tamM[k - 1] + 1):
                    if rl.relM[i][j][k][l][0][0] != 0:
                        res = rl.relM[i][j][k][l]
                        cad = "[Z" + str(i) + str(j) + ", Z" + str(k) + str(l) + "] - ("
                        primer = True
                        resq = sp.S(0)
                        for m in res:
                            resq += sp.Rational(m[0], m[1])*Z[m[2]][m[3]]
                            if m[0] > 0 and primer == False:
                                cad += "+ "
                            elif m[0] < 0:
                                cad += "- "
                            primer = False
                            if abs(m[0]) > 1:
                                cad += str(abs(m[0])) + "*"
                            cad += "Z" + str(m[2]) + str(m[3]) + ""
                            if m[1] != 1:
                                cad += "/" + str(m[1])
                            cad += " "
                        cad = cad[:-1] + ") ="
                        r = (Commutator(Z[i][j], Z[k][l]) - resq).doit().expand()
                        if (r != 0):
                            print(cad, r)
