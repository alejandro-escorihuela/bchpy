#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 21-12-2020
# alex
# vmod.py

import numpy as np
import sympy as sp
import time as tm
import sys
sys.path.insert(0, '../')
from bchpy import *
import relations as rl

if __name__ == "__main__":
    ord_bch = 4
    c = sp.Symbol("c")
    a = sp.symbols("a0:4", real=True)
    b = sp.symbols("b0:4")
    # print("Vmod 1")
    # t0 = tm.time()
    # print("1a iteració")
    # esq = bch9(a[0]*Eel(1, 1), b[0]*Eel(1, 2), depth = ord_bch)
    # print("2a iteració")
    # esq = bch9(b[1]*(Eel(1, 2) + c*Eel(3, 2)), esq, depth = ord_bch)
    # print("3a iteració")
    # esq = bch9(a[0]*Eel(1, 1), esq, depth = ord_bch)
    # print("4a iteració")
    # esq = bch9(b[0]*Eel(1, 2), esq, depth = ord_bch)
    # print("Conversió a mètode")
    # met = Metode(ord_bch)
    # met.importFromExpr(esq)
    # met.cprint()
    # t1 = tm.time()
    # print(t1 - t0, "s")
    # A = sp.nsolve((met.w[1][1] - 1, met.w[1][2] - 1, met.w[3][1], met.w[3][2]), (a[0], b[0], b[1], c), (1, 1, 1 ,1))
    # print(A)

    # print("Vmod 2")
    # t0 = tm.time()
    # print("1a iteració")
    # esq = bch9(a[0]*Eel(1, 1), b[0].conjugate().expand(complex=True)*Eel(1, 2), depth = ord_bch)
    # print("2a iteració")
    # esq = bch9(b[1]*(Eel(1, 2) + c*Eel(3, 2)), esq, depth = ord_bch)
    # print("3a iteració")
    # esq = bch9(a[0]*Eel(1, 1), esq, depth = ord_bch)
    # print("4a iteració")
    # esq = bch9(b[0].expand(complex=True)*Eel(1, 2), esq, depth = ord_bch)
    # print("Conversió a mètode")
    # met = Metode(ord_bch)
    # met.importFromExpr(esq)
    # met.cprint()
    # t1 = tm.time()
    # print(t1 - t0, "s")
    # A = sp.nsolve((met.w[1][1] - 1, met.w[1][2] - 1, met.w[2][1], met.w[3][1], met.w[3][2]), (a[0], re(b[0]), im(b[0]), b[1], c), (1, 1, 1, 1, 1))
    # print(A)

    # print("Vmod 3")
    # t0 = tm.time()
    # print("1a iteració")
    # esq = bch9(b[0].conjugate().expand(complex=True)*(Eel(1, 2) + c*Eel(3, 2)), a[0]*Eel(1,1), depth = ord_bch)
    # print("2a iteració")
    # esq = bch9(a[1]*Eel(1, 1), esq, depth = ord_bch)
    # print("3a iteració")
    # esq = bch9(b[0].expand(complex=True)*(Eel(1, 2) + c*Eel(3, 2)), esq, depth = ord_bch)
    # print("4a iteració")
    # esq = bch9(a[0]*Eel(1,1), esq, depth = ord_bch)
    # print("Conversió a mètode")
    # met = Metode(ord_bch)
    # met.importFromExpr(esq)
    # met.cprint()
    # t1 = tm.time()
    # print(t1 - t0, "s")
    # A = sp.nsolve((met.w[1][1] - 1, met.w[1][2] - 1, met.w[2][1], met.w[3][1], met.w[3][2]), (a[0], re(b[0]), im(b[0]), b[1], c), (1, 1, 1, 1, 1))
    # print(A)
    
