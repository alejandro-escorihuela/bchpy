#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 09-12-2020
# alex
# recursives.py

import numpy as np
import sympy as sp
import time as tm
import sys
sys.path.insert(0, '../')
from bchpy import *
import relations as rl

if __name__ == "__main__":
    ord_bch = 8
    print("exp(bet(i,j)*Eij)=exp(x*E12)*exp(alp(i,j)*Eij)")
    x = sp.Symbol("x")
    t0 = tm.time()
    esq = sp.S(0)
    alp = sp.MatrixSymbol("alp", ord_bch + 1, rl.tamE[ord_bch - 1] + 1)
    for i in range(ord_bch):
        for j in range(rl.tamE[i]):
            esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1)
    #esq = bch9(x*Eel(1, 2), esq, depth = ord_bch)
    esq = bch9(Eel(1, 1), Eel(1, 2), depth = ord_bch, debug = True)
    esq = bch9(Eel(1, 2), esq, depth = ord_bch, debug = True)
    metBD = Metode(ord_bch)
    metBD.importFromExpr(esq, debug = True)
    metBD.save("recB.dat")
    print(metBD)
    t1 = tm.time()
    print(t1 - t0, "s")
    print("\n")
    
    # print("exp(gam(i, j)*Eij)=exp(y*E11)*exp(alp(i, j)*Eij)")
    # y = sp.Symbol("y")
    # t0 = tm.time()
    # esq = sp.S(0)
    # alp = sp.MatrixSymbol("alp", ord_bch + 1, rl.tamE[ord_bch - 1] + 1)
    # for i in range(ord_bch):
    #     for j in range(rl.tamE[i]):
    #         esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1)
    # esq = bch9(y*Eel(1, 1), esq, depth = ord_bch)
    # metAD = Metode(ord_bch)
    # metAD.importFromExpr(esq)
    # print(metAD)
    # t1 = tm.time()
    # print(t1 - t0, "s")
