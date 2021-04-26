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
    ord_bch = 5
    x, y = sp.symbols("x y")
    t0 = tm.time()
    esq = sp.S(0)
    xa = sp.S(0)
    xb = sp.S(0)
    alp = sp.MatrixSymbol("alp", ord_bch + 1, rl.tamZ[ord_bch - 1] + 1)
    for i in range(ord_bch):
        xa += x**(i + 1)*Eel(i + 1, 1, "Z")
        xb += (-1)**i*y**(i + 1)*Eel(i + 1, 1, "Z")
        for j in range(rl.tamZ[i]):
            esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1, "Z")
    print(xa)
    print(xb)
    #esq = bch9(xa, esq, depth = ord_bch, debug = True)
    esq = bch9(xa, xb, depth = ord_bch, debug = True)
    print(esq)
    metBD = Metode(ord_bch, "Z")
    metBD.importFromExpr(esq, debug = True)
    print(metBD)
    t1 = tm.time()
    print(t1 - t0, "s")
    
    # print("exp(w(i,j)*Eij)=exp(x*E12)*exp(alp(i,j)*Eij)")
    # x, y = sp.symbols("x y")
    # t0 = tm.time()
    # esq = sp.S(0)
    # alp = sp.MatrixSymbol("alp", ord_bch + 1, rl.tamE[ord_bch - 1] + 1)
    # for i in range(ord_bch):
    #     for j in range(rl.tamE[i]):
    #         esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1)
    # esq = bch9(x*Eel(1, 2), esq, depth = ord_bch, debug = True)
    # metBD = Metode(ord_bch)
    # metBD.importFromExpr(esq, debug = True)
    # print(metBD)
    # t1 = tm.time()
    # print(t1 - t0, "s")
    
    # print("exp(w(i, j)*Eij)=exp(y*E11)*exp(alp(i, j)*Eij)")
    # x = sp.Symbol("x")
    # t0 = tm.time()
    # esq = sp.S(0)
    # alp = sp.MatrixSymbol("alp", ord_bch + 1, rl.tamE[ord_bch - 1] + 1)
    # for i in range(ord_bch):
    #     for j in range(rl.tamE[i]):
    #         esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1)
    # esq = bch9(x*Eel(1, 1), esq, depth = ord_bch, debug = True)
    # metAD = Metode(ord_bch)
    # metAD.importFromExpr(esq, debug = True)
    # print(metAD)
    # t1 = tm.time()
    # print(t1 - t0, "s")

    # print("exp(w(i,j)*Eij)=exp(x*E12)*exp(alp(i,j)*Eij)*exp(x*E12)")
    # x, y = sp.symbols("x y")
    # t0 = tm.time()
    # esq = sp.S(0)
    # alp = sp.MatrixSymbol("alp", ord_bch + 1, rl.tamE[ord_bch - 1] + 1)
    # for i in range(0, ord_bch, 2):
    #     for j in range(rl.tamE[i]):
    #         esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1)
    # esq = bch9(esq, x*Eel(1, 2), depth = ord_bch, debug = True)
    # esq = bch9(x*Eel(1, 2), esq.expand(), depth = ord_bch, debug = True)
    # metBDB = Metode(ord_bch)
    # metBDB.importFromExpr(esq, debug = True)
    # print(metBDB)
    # t1 = tm.time()
    # print(t1 - t0, "s")  
    
    # print("exp(w(i,j)*Eij)=exp(x*E11)*exp(alp(i,j)*Eij)*exp(x*E11)")
    # x, y = sp.symbols("x y")
    # t0 = tm.time()
    # esq = sp.S(0)
    # alp = sp.MatrixSymbol("alp", ord_bch + 1, rl.tamE[ord_bch - 1] + 1)
    # for i in range(0, ord_bch, 2):
    #     for j in range(rl.tamE[i]):
    #         esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1)
    # esq = bch9(esq, x*Eel(1, 1), depth = ord_bch, debug = True)
    # esq = bch9(x*Eel(1, 1), esq.expand(), depth = ord_bch, debug = True)
    # metBDB = Metode(ord_bch)
    # metBDB.importFromExpr(esq, debug = True)
    # print(metBDB)
    # t1 = tm.time()
    # print(t1 - t0, "s") 

    # print("exp(w(i,j)*Eij)=exp(x*E12)*exp(alp(i,j)*Eij)*exp(x.conjugate()*E12)")
    # x, y = sp.symbols("x y")
    # t0 = tm.time()
    # esq = sp.S(0)
    # alp = sp.MatrixSymbol("alp", ord_bch + 1, rl.tamE[ord_bch - 1] + 1)
    # for i in range(ord_bch):
    #     for j in range(rl.tamE[i]):
    #         if (i + 1) % 2 == 0:
    #             esq += sp.I*alp[i + 1, j + 1]*Eel(i + 1, j + 1)
    #         else:
    #             esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1)
    # esq = bch9(esq, x.conjugate().expand(complex=True)*Eel(1, 2), depth = ord_bch, debug = True)
    # esq = collectEsq(esq)
    # esq = bch9(x.expand(complex=True)*Eel(1, 2), esq, depth = ord_bch, debug = True)
    # metBCB = Metode(ord_bch)
    # metBCB.importFromExpr(esq, debug = True)
    # metBCB.collectI()
    # print(metBCB)
    # t1 = tm.time()
    # print(t1 - t0, "s")  
    

    # print("exp(w(i,j)*Eij)=exp(x*E11)*exp(alp(i,j)*Eij)*exp(x.conjugate()*E11)")
    # x, y = sp.symbols("x y")
    # t0 = tm.time()
    # esq = sp.S(0)
    # alp = sp.MatrixSymbol("alp", ord_bch + 1, rl.tamE[ord_bch - 1] + 1)
    # for i in range(ord_bch):
    #     for j in range(rl.tamE[i]):
    #         if (i + 1) % 2 == 0:
    #             esq += sp.I*alp[i + 1, j + 1]*Eel(i + 1, j + 1)
    #         else:
    #             esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1)
    # esq = bch9(esq, x.conjugate().expand(complex=True)*Eel(1, 1), depth = ord_bch, debug = True)
    # esq = collectEsq(esq)
    # esq = bch9(x.expand(complex=True)*Eel(1, 1), esq, depth = ord_bch, debug = True)
    # metACA = Metode(ord_bch)
    # metACA.importFromExpr(esq, debug = True)
    # metACA.collectI()
    # print(metACA)
    # t1 = tm.time()
    # print(t1 - t0, "s")       
