#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 09-12-2020
# alex
# recursives.py

import numpy as np
import sympy as sp
import relations as rl
import time as tm
from bchpy import *

if __name__ == "__main__":
    print("exp(bet(i,j)*Eij)=exp(x*E12)*exp(alp(i,j)*Eij)")
    x = sp.Symbol("x")
    t0 = tm.time()
    esq = sp.S(0)
    alp = sp.MatrixSymbol("alp", len(rl.tamE) + 1, rl.tamE[-1] + 1)
    for i in range(len(rl.tamE)):
        for j in range(rl.tamE[i]):
            esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1)
    esq = bch6(x*Eel(1,2), esq)
    metBD = Metode()
    metBD.importFromExpr(esq)
    print(metBD)
    t1 = tm.time()
    print(t1 - t0, "s")
    print("\n")
    
    print("exp(gam(i,j)*Eij)=exp(y*E11)*exp(alp(i,j)*Eij)")
    y = sp.Symbol("y")
    t0 = tm.time()
    esq = sp.S(0)
    for i in range(len(rl.tamE)):
        for j in range(rl.tamE[i]):
            esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1)
    esq = bch6(y*Eel(1,1), esq)
    metAD = Metode()
    metAD.importFromExpr(esq)
    print(metAD)
    t1 = tm.time()
    print(t1 - t0, "s")
