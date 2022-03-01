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

def printmetC(met, o):
    cad = []
    v = rl.tamE
    tam = (v[o - 1] + 1)
    for i in range(len(met.w)):
        for j in range(len(met.w[i])):
            cad.append(str(sp.ccode(met.w[i][j])))
    for k in range(len(cad)):
        for t in range(tam, tam**2 + tam):
            p = t - tam
            j = p%tam
            i = p//tam + 1
            l = sum(v[0:i - 1:2]) + j - 1
            cad[k] = cad[k].replace("alp[" + str(t) + "]", "alp[" + str(l) + "]")
    oset = 2
    ind = 0
    for k in range(oset + 1, oset + 1 + sum(v[:1])):
        print("res[%d] = %s;" % (ind, cad[k]))
        ind += 1    
    for l in range(3, o + 1, 2):
        for k in range(oset + 1 + (l - 1) + sum(v[:l - 1]), oset + l + sum(v[:l])):
            print("res[%d] = %s;" % (ind, cad[k]))
            ind += 1
            
if __name__ == "__main__":
    ord_bch = 9
    # print("?¿?")
    # x, y = sp.symbols("x y")
    # t0 = tm.time()
    # esq = sp.S(0)
    # xa = sp.S(0)
    # xb = sp.S(0)
    # alp = sp.MatrixSymbol("alp", ord_bch + 1, rl.tamZ[ord_bch - 1] + 1)
    # for i in range(ord_bch):
    #     xa += x**(i + 1)*Eel(i + 1, 1, "Z")
    #     xb += (-1)**i*x**(i + 1)*Eel(i + 1, 1, "Z")
    #     for j in range(rl.tamZ[i]):
    #         esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1, "Z")

    # esq = bch9(xb, esq, depth = ord_bch, debug = True)
    # metBD = Metode(ord_bch, "Z")
    # metBD.importFromExpr(esq, debug = True)
    # print(metBD)
    # t1 = tm.time()
    # print(t1 - t0, "s")
    
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
    #         esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1, "E")
    # esq = bch9(esq, x*Eel(1, 2, "E"), depth = ord_bch, debug = True)
    # esq = bch9(x*Eel(1, 2, "E"), esq.expand(), depth = ord_bch, debug = True)
    # metBDB = Metode(ord_bch, "E")
    # metBDB.importFromExpr(esq, debug = True)
    # printmetC(metBDB, ord_bch)
    # t1 = tm.time()
    # print(t1 - t0, "s")  
    
    print("exp(w(i,j)*Eij)=exp(x*E11)*exp(alp(i,j)*Eij)*exp(x*E11)")
    x, y = sp.symbols("x y")
    t0 = tm.time()
    esq = sp.S(0)
    alp = sp.MatrixSymbol("alp", ord_bch + 1, rl.tamE[ord_bch - 1] + 1)
    for i in range(0, ord_bch, 2):
        for j in range(rl.tamE[i]):
            esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1, "E")
    esq = bch9(esq, x*Eel(1, 1, "E"), depth = ord_bch, debug = True)
    esq = bch9(x*Eel(1, 1, "E"), esq.expand(), depth = ord_bch, debug = True)
    metADA = Metode(ord_bch)
    metADA.importFromExpr(esq, debug = True)
    #print(metADA)
    printmetC(metADA, ord_bch)
    t1 = tm.time()
    print(t1 - t0, "s") 

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
