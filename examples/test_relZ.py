#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 24-04-2021
# alex
# test_relZ.py

import sympy as sp
from sympy.physics.quantum import Operator, Commutator
import sys
sys.path.insert(0, '../')
import relations as rl
from qolib import *

if __name__ == "__main__":
    n_max = 9
    Y = []
    Y.append([0])
    for i in range(n_max + 1):
        Y.append(Operator("Y_{%d}" % (i + 1)))
        
    Z = llegir_baseCAdj("compAdj_base.txt", Y)
    
    ord_bch = len(rl.tamZ) + 1
    for i in range(1, ord_bch):
        for j in range(1, rl.tamZ[i - 1] + 1):
            for k in range(1, ord_bch):
                for l in range(1, rl.tamZ[k - 1] + 1):
                    if rl.relZ[i][j][k][l][0][0] != 0:
                        res = rl.relZ[i][j][k][l]
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
