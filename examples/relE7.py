#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 28-12-2020
# alex
# relE7.py

import sympy as sp
from sympy.physics.quantum import Operator, Commutator
import itertools
import sys
sys.path.insert(0, '../')
from bchpy import *
import relations as rl

def resoldre(exprEsq, exprDre):
    eqs = []
    for i in exprDre.args:
        cnc_i = i.args_cnc()
        equ = sp.S(0)
        for j in exprEsq.args:
            cnc_j = j.args_cnc()
            if cnc_i[1] == cnc_j[1]:
                equ_esq = sp.S(1)
                for elem in cnc_j[0]:
                    equ_esq *= elem
                equ += equ_esq
        equ_dre = sp.S(1)
        for elem in cnc_i[0]:
            equ_dre *= elem
        equ -= equ_dre
        if equ not in eqs:
            eqs.append(equ)
    return sp.solve(eqs)

            
if __name__ == "__main__":
    A = Operator("A")
    B = Operator("B")
    c = sp.symbols("c1:19")
    E = [[0 for j in range(19)] for i in range(8)]
    E[1][1] = A
    E[1][2] = B
    E[2][1] = Commutator(A,B)
    E[3][1] = Commutator(A, Commutator(A, B))
    E[3][2] = Commutator(B, Commutator(A, B))
    E[4][1] = Commutator(A, Commutator(A, Commutator(A, B)))
    E[4][2] = Commutator(B, Commutator(A, Commutator(A, B)))
    E[4][3] = Commutator(B, Commutator(B, Commutator(B, A)))
    E[5][1] = Commutator(A, Commutator(A, Commutator(A, Commutator(A, B))))
    E[5][2] = Commutator(B, Commutator(A, Commutator(A, Commutator(A, B))))
    E[5][3] = Commutator(A, Commutator(A, Commutator(B, Commutator(B, A))))
    E[5][4] = Commutator(B, Commutator(B, Commutator(A, Commutator(A, B))))
    E[5][5] = Commutator(A, Commutator(B, Commutator(B, Commutator(B, A))))
    E[5][6] = Commutator(B, Commutator(B, Commutator(B, Commutator(B, A))))
    E[6][1] = Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B)))))
    E[6][2] = Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B)))))
    E[6][3] = Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B)))))
    E[6][4] = Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, B)))))
    E[6][5] = Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B)))))
    E[6][6] = Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(B, A)))))
    E[6][7] = Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(B, A)))))
    E[6][8] = Commutator(A, Commutator(B, Commutator(B, Commutator(B, Commutator(B, A)))))
    E[6][9] = Commutator(B, Commutator(B, Commutator(B, Commutator(B, Commutator(B, A)))))
    for i in range(1, 10):
        E[7][2*i - 1] = Commutator(A, E[6][i])  
        E[7][2*i] = Commutator(B, E[6][i])
    cesq = sp.S(0)
    for i in range(1, 19):
        cesq += c[i - 1]*E[7][i]
    cesq = cesq.doit().expand()
    for ind in ([[1, 6], [2, 5], [3, 4]]):
        i, k = ind
        for j, l in itertools.product(range(1, rl.tamE[i - 1] + 1), range(1, rl.tamE[k - 1] + 1)):
            cdre = Commutator(E[i][j], E[k][l]).doit().expand()
            sol = resoldre(cesq, cdre)
            cad = "[E" + str(i) + str(j) + ", E" + str(k) + str(l) + "] = "
            cadrel = "relE[" + str(i) + "][" + str(j) + "][" + str(k) + "][" + str(l) + "] = ["
            for m in sol:
                if sol[m] != 0:
                    if sol[m] < 0:
                        cad += str(sol[m]) + "*E7" + str(c.index(m) + 1)
                    else:
                        cad += "+" + str(sol[m]) + "*E7" + str(c.index(m) + 1)
                    nume, deno = sp.fraction(sol[m])
                    cadrel += "[" + str(nume) + ", " + str(deno) + ", 7, " + str(c.index(m) + 1) + "], "
            cadrel = cadrel[:-2] + "]"
            # print(cad)
            print(cadrel)
