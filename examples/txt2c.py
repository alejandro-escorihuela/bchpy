#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 08-02-2024
# alex
# txt2c.py

import sys
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr
sys.path.insert(0, '../')
from bchpy import *
import relations as rl

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Hi ha que passar almenys un fitxer amb les solucions.")
        exit(-1)
    ord_bch = 9
    x, y = sp.symbols("x y")
    esq = sp.S(0)
    xa = sp.S(0)
    xb = sp.S(0)
    alp = sp.MatrixSymbol("alp", ord_bch + 1, rl.tamC[ord_bch - 1] + 1)
    for i in range(ord_bch):
        xa += x**(i + 1)*Eel(i + 1, 1, "C")
        xb += (-1)**i*x**(i + 1)*Eel(i + 1, 1, "C")
        for j in range(rl.tamC[i]):
            esq += alp[i + 1, j + 1]*Eel(i + 1, j + 1, "C")    
    nom_fit = sys.argv[1]
    f = open(nom_fit, "r")
    linies = f.readlines()
    cad = []
    for lin in linies:
        cad.append(sp.ccode(eval(lin.split(" = ")[1].replace("\n", ""))))
    tam = (rl.tamC[ord_bch - 1] + 1)    
    for k in range(len(cad)):
        for t in range(tam, tam**2 + tam):
            # t, p, j, i, l = valor vell, valor vell sense tenir en compte la fila de zeros, columna, fila, valor nou
            p = t - tam
            j = p%tam
            i = p//tam + 1
            l = sum(rl.tamC[0:i - 1:1]) + j - 1
            # print("t p j i l = ", t, p, j, i , l)
            cad[k] = cad[k].replace("alp[" + str(t) + "]", "alp[" + str(l) + "]")    
    for k in range(len(cad)):
        print("res[%d] = %s;" % (k, cad[k]))
            
