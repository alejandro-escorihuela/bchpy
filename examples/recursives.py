#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 09-12-2020
# alex
# recursives.py
# per provar BCH

import numpy as np
import sympy as sp
import time as tm
import sys
sys.path.insert(0, '../')
from bchpy import *
import relations as rl

def printmetC(met, o, simconj = False, cq = 2, base = "E"):
    cad = []

    if base == "E":
        v = rl.tamE
    elif base == "C":
        v = rl.tamC
    elif base == "M":
        v = rl.tamM
    elif base == "M4":
        v = rl.tamM4
    elif base == "M6":
        v = rl.tamM6        
    tam = (v[o - 1] + 1)
    for i in range(len(met.w)):
        for j in range(len(met.w[i])):
            if simconj:
                met.w[i][j] = met.w[i][j].subs(re(x), xre).subs(im(x), xim).subs(sp.I, sp.S(1))
            cad.append(str(sp.ccode(met.w[i][j])))
    # si cq = 1 es tenen en compte els elements d'ordre parell
    # si cq = 2 no es tenen en compte
    lim = tam**2 + tam
    if (base == "M" or base == "M4") and lim <= 100:
        lim = 100
    elif base == "M6" and lim <= 1000:
        lim = 1000   
    for k in range(len(cad)):
        for t in range(tam, lim):
            # print(k, "->", cad[k])
            # t, p, j, i, l = valor vell, valor vell sense tenir en compte la fila de zeros, columna, fila, valor nou
            p = t - tam
            j = p%tam
            i = p//tam + 1
            l = sum(v[0:i - 1:cq]) + j - 1
            # print("t p j i l = ", t, p, j, i , l)
            cad[k] = cad[k].replace("alp[" + str(t) + "]", "alp[" + str(l) + "]")
    ind = 0
    if not simconj:
        oset = 2
        for k in range(oset + 1, oset + 1 + sum(v[:1])):
            print("res[%d] = %s;" % (ind, cad[k]))
            ind += 1
        for l in range(1 + cq, o + 1, cq):
            for k in range(oset + 1 + (l - 1) + sum(v[:l - 1]), oset + l + sum(v[:l])):
                print("res[%d] = %s;" % (ind, cad[k]))
                ind += 1
    else:
        for k in range(len(cad)):
            if cad[k] != "0":
                print("res[%d] = %s;" % (ind, cad[k]))
                ind += 1

if __name__ == "__main__":
    ord_bch = 9
    
    # print("Mètode i adjunt")
    # x, y = sp.symbols("x y")
    # t0 = tm.time()
    # esq = sp.S(0)
    # xa = sp.S(0)
    # xb = sp.S(0)
    # alp = sp.MatrixSymbol("alp", ord_bch + 1, rl.tamC[ord_bch - 1] + 1)
    # for i in range(ord_bch):
    #     xa += x**(i + 1)*Eel(i + 1, 1, "C")
    #     xb += (-1)**i*x**(i + 1)*Eel(i + 1, 1, "C")
    #     for j in range(rl.tamC[i]):
    #         esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1, "C")

    # esq = bch20(xa, esq, depth = ord_bch, debug = True)
    # metBD = Metode(ord_bch, "C")
    # metBD.importFromExpr(esq, debug = True)
    # # print(xa)
    # # print(xb)
    # # print(metBD)
    # printmetC(metBD, ord_bch, simconj = False, cq = 1, base = "C")
    # t1 = tm.time()
    # print(t1 - t0, "s")

    print("Mètode i adjunt simètric")
    x, y = sp.symbols("x y")
    t0 = tm.time()
    esq = sp.S(0)
    xa = sp.S(0)
    xb = sp.S(0)
    alp = sp.MatrixSymbol("alp", ord_bch + 1, rl.tamC[ord_bch - 1] + 1)
    for i in range(ord_bch):
        xa += x**(i + 1)*Eel(i + 1, 1, "C")
        xb += (-1)**i*x**(i + 1)*Eel(i + 1, 1, "C")
        for j in range(rl.tamC[i]):
            esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1, "C")

    # esq = bch20(esq, xa, depth = ord_bch, debug = True)
    # esq = bch20(xb, esq.expand(), depth = ord_bch, debug = True)
    # metBD = Metode(ord_bch, "C")
    # metBD.importFromExpr(esq, debug = True)
    # # print(xa)
    # # print(xb)
    # # print(metBD)
    # printmetC(metBD, ord_bch, simconj = False, cq = 1, base = "C")
    # t1 = tm.time()
    # print(t1 - t0, "s")
    
    # print("exp(w(i,j)*Mij)=exp(x*M11+x**3*M31+x**5*M51+...)*exp(alp(i,j)*Mij)*exp(x*M11+x**3*M31+x**5*M51+...)")
    # x, y = sp.symbols("x y")
    # t0 = tm.time()
    # esq = sp.S(0)
    # xa = sp.S(0)
    # xb = sp.S(0)
    # ss = x*Eel(1, 1, "M")
    # for i in range(2, ord_bch, 2):
    #     ss +=x**(i + 1)*Eel(i + 1, 1, "M")
    # alp = sp.MatrixSymbol("alp", ord_bch + 1, rl.tamM[ord_bch - 1] + 1)
    # for i in range(0, ord_bch, 2):
    #     for j in range(rl.tamM[i]):
    #         esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1, "M")
    # print("Z   =", esq)
    # print("S*x =", ss)
    # esq = bch20(esq, ss, depth = ord_bch, debug = True)
    # esq = bch20(ss, esq.expand(), depth = ord_bch, debug = True)
    # metBD = Metode(ord_bch, "M")
    # metBD.importFromExpr(esq, debug = True)
    # # print(metBD)
    # printmetC(metBD, ord_bch, simconj = False, cq = 2, base = "M")
    # t1 = tm.time()
    # print(t1 - t0, "s")

    # print("exp(w(i,j)*Mij)=exp(alp(i,j)*Mij)*exp(x*M11+x**3*M31+x**5*M51+...)")
    # x, y = sp.symbols("x y")
    # t0 = tm.time()
    # esq = sp.S(0)
    # ss = x*Eel(1, 1, "M")
    # for i in range(2, ord_bch, 2):
    #     ss +=x**(i + 1)*Eel(i + 1, 1, "M")
    # alp = sp.MatrixSymbol("alp", ord_bch + 1, rl.tamM[ord_bch - 1] + 1)
    # for i in range(ord_bch):
    #     for j in range(rl.tamM[i]):
    #         esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1, "M")
    # print("Z   =", esq)
    # print("S*x =", ss)
    # esq = bch20(ss, esq, depth = ord_bch, debug = True)
    # metBD = Metode(ord_bch, "M")
    # metBD.importFromExpr(esq, debug = True)
    # # print(metBD)
    # printmetC(metBD, ord_bch, simconj = False, cq = 1, base = "M")
    # t1 = tm.time()
    # print(t1 - t0, "s")
    
    # print("exp(w(i,j)*Mij)=exp(x*M11+x**5*M51+x**7*M71+...)*exp(alp(i,j)*Mij)*exp(x*M11+x**5*M51+x**7*M71+...)")
    # x, y = sp.symbols("x y")
    # t0 = tm.time()
    # esq = sp.S(0)
    # xa = sp.S(0)
    # xb = sp.S(0)
    # ss = x*Eel(1, 1, "M4") + x**5*Eel(5, 1, "M4")
    # for i in range(6, ord_bch, 2):
    #     ss +=x**(i + 1)*Eel(i + 1, 1, "M4")
    # alp = sp.MatrixSymbol("alp", ord_bch + 1, rl.tamM4[ord_bch - 1] + 1)
    # for i in range(0, ord_bch, 2):
    #     for j in range(rl.tamM4[i]):
    #         esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1, "M4")
    # print("Z   =", esq)
    # print("S*x =", ss)
    # esq = bch20(esq, ss, depth = ord_bch, debug = True)
    # esq = bch20(ss, esq.expand(), depth = ord_bch, debug = True)
    # metBD = Metode(ord_bch, "M4")
    # metBD.importFromExpr(esq, debug = True)
    # print(metBD)
    # printmetC(metBD, ord_bch, simconj = False, cq = 2, base = "M4")
    # t1 = tm.time()
    # print(t1 - t0, "s")
    
    # print("exp(w(i,j)*Mij)=exp(alp(i,j)*Mij)*exp(x*M11+x**5*M51+x**7*M71+...)")
    # x, y = sp.symbols("x y")
    # t0 = tm.time()
    # esq = sp.S(0)
    # ss = x*Eel(1, 1, "M4") + x**5*Eel(5, 1, "M4")
    # for i in range(6, ord_bch, 2):
    #     ss +=x**(i + 1)*Eel(i + 1, 1, "M4")
    # alp = sp.MatrixSymbol("alp", ord_bch + 1, rl.tamM4[ord_bch - 1] + 1)
    # for i in range(ord_bch):
    #     for j in range(rl.tamM4[i]):
    #         esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1, "M4")
    # print("Z       =", esq)
    # print("S*x     =", ss) 
    # esq = bch20(ss, esq, depth = ord_bch, debug = True)
    # metBD = Metode(ord_bch, "M4")
    # metBD.importFromExpr(esq, debug = True)
    # print(metBD)
    # printmetC(metBD, ord_bch, simconj = False, cq = 1, base = "M4")
    # t1 = tm.time()
    # print(t1 - t0, "s")
    
    # print("exp(w(i,j)*Mij)=exp(x*M11+x**7*M71+x**9*M91+...)*exp(alp(i,j)*Mij)*exp(x*M11+x**7*M71+x**9*M91+...)")
    # x, y = sp.symbols("x y")
    # t0 = tm.time()
    # esq = sp.S(0)
    # xa = sp.S(0)
    # xb = sp.S(0)
    # ss = x*Eel(1, 1, "M6") + x**7*Eel(7, 1, "M6")
    # for i in range(8, ord_bch, 2):
    #     ss +=x**(i + 1)*Eel(i + 1, 1, "M6")
    # alp = sp.MatrixSymbol("alp", ord_bch + 1, rl.tamM6[ord_bch - 1] + 1)
    # for i in range(0, ord_bch, 2):
    #     for j in range(rl.tamM6[i]):
    #         esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1, "M6")
    # print("Z   =", esq)
    # print("S*x =", ss)
    # esq = bch20(esq, ss, depth = ord_bch, debug = True)
    # esq = bch20(ss, esq.expand(), depth = ord_bch, debug = True)
    # metBD = Metode(ord_bch, "M6")
    # metBD.importFromExpr(esq, debug = True)
    # print(metBD)
    # printmetC(metBD, ord_bch, simconj = False, cq = 2, base = "M6")
    # t1 = tm.time()
    # print(t1 - t0, "s")

    # print("exp(w(i,j)*Mij)=exp(alp(i,j)*Mij)*exp(x*M11+x**7*M71+x**9*M91+...)")
    # x, y = sp.symbols("x y")
    # t0 = tm.time()
    # esq = sp.S(0)
    # xa = sp.S(0)
    # xb = sp.S(0)
    # ss = x*Eel(1, 1, "M6") + x**7*Eel(7, 1, "M6")
    # for i in range(8, ord_bch, 2):
    #     ss +=x**(i + 1)*Eel(i + 1, 1, "M6")
    # alp = sp.MatrixSymbol("alp", ord_bch + 1, rl.tamM6[ord_bch - 1] + 1)
    # for i in range(0, ord_bch, 2):
    #     for j in range(rl.tamM6[i]):
    #         esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1, "M6")
    # print("Z   =", esq)
    # print("S*x =", ss)
    # esq = bch20(ss, esq, depth = ord_bch, debug = True)
    # metBD = Metode(ord_bch, "M6")
    # metBD.importFromExpr(esq, debug = True)
    # print(metBD)
    # printmetC(metBD, ord_bch, simconj = False, cq = 1, base = "M6")
    # t1 = tm.time()
    # print(t1 - t0, "s")
    
    # print("exp(w(i,j)*Eij)=exp(x*E12)*exp(alp(i,j)*Eij)")
    # x, y = sp.symbols("x y")
    # t0 = tm.time()
    # esq = sp.S(0)
    # alp = sp.MatrixSymbol("alp", ord_bch + 1, rl.tamE[ord_bch - 1] + 1)
    # for i in range(ord_bch):
    #     for j in range(rl.tamE[i]):
    #         esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1, "E")
    # esq = bch20(x*Eel(1, 2, "E"), esq, depth = ord_bch, debug = True)
    # metBD = Metode(ord_bch)
    # metBD.importFromExpr(esq, debug = True)
    # printmetC(metBD, ord_bch, simconj = False, cq = 1)
    # # print(metBD)
    # t1 = tm.time()
    # print(t1 - t0, "s")
    
    # print("exp(w(i, j)*Eij)=exp(x*E11)*exp(alp(i, j)*Eij)")
    # x = sp.Symbol("x")
    # t0 = tm.time()
    # esq = sp.S(0)
    # alp = sp.MatrixSymbol("alp", ord_bch + 1, rl.tamE[ord_bch - 1] + 1)
    # for i in range(ord_bch):
    #     for j in range(rl.tamE[i]):
    #         esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1, "E")
    # esq = bch20(x*Eel(1, 1, "E"), esq, depth = ord_bch, debug = True)
    # metAD = Metode(ord_bch)
    # metAD.importFromExpr(esq, debug = True)
    # printmetC(metAD, ord_bch, simconj = False, cq = 1)
    # # print(metAD)
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
    
    # print("exp(w(i,j)*Eij)=exp(x*E11)*exp(alp(i,j)*Eij)*exp(x*E11)")
    # x, y = sp.symbols("x y")
    # t0 = tm.time()
    # esq = sp.S(0)
    # alp = sp.MatrixSymbol("alp", ord_bch + 1, rl.tamE[ord_bch - 1] + 1)
    # for i in range(0, ord_bch, 2):
    #     for j in range(rl.tamE[i]):
    #         esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1, "E")
    # esq = bch20(esq, x*Eel(1, 1, "E"), depth = ord_bch, debug = True)
    # esq = bch20(x*Eel(1, 1, "E"), esq.expand(), depth = ord_bch, debug = True)
    # metADA = Metode(ord_bch)
    # metADA.importFromExpr(esq, debug = True)
    # #print(metADA)
    # printmetC(metADA, ord_bch)
    # t1 = tm.time()
    # print(t1 - t0, "s") 

