#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 28-12-2020
# alex
# niats.py
# Corxets niats o niuats

import sympy as sp
from sympy.physics.quantum import Operator, Commutator
import itertools
import sys
sys.path.insert(0, '../')
from qolib import *

if __name__ == "__main__":
    n_max = 10
    A = Operator("A")
    B = Operator("B")

    E = []
    E.append([0])

    E.append([0])
    E[1].append(A)
    E[1].append(B)

    E.append([0])
    E[2].append(Commutator(A, B))
    
    E.append([0])
    E[3].append(Commutator(A, E[2][1]))
    E[3].append(Commutator(B, E[2][1]))

    E.append([0])
    E[4].append(Commutator(A, Commutator(A, Commutator(A, B))))
    E[4].append(Commutator(B, Commutator(A, Commutator(A, B))))
    E[4].append(Commutator(B, Commutator(B, Commutator(B, A))))

    E.append([0])
    E[5].append(Commutator(A, Commutator(A, Commutator(A, Commutator(A, B)))))
    E[5].append(Commutator(B, Commutator(A, Commutator(A, Commutator(A, B)))))
    E[5].append(Commutator(A, Commutator(A, Commutator(B, Commutator(B, A)))))
    E[5].append(Commutator(B, Commutator(B, Commutator(A, Commutator(A, B)))))
    E[5].append(Commutator(A, Commutator(B, Commutator(B, Commutator(B, A)))))
    E[5].append(Commutator(B, Commutator(B, Commutator(B, Commutator(B, A)))))

    E.append([0])
    E[6].append(Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B))))))
    E[6].append(Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B))))))
    E[6].append(Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B))))))
    E[6].append(Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, B))))))
    E[6].append(Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B))))))
    E[6].append(Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(B, A))))))
    E[6].append(Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(B, A))))))
    E[6].append(Commutator(A, Commutator(B, Commutator(B, Commutator(B, Commutator(B, A))))))
    E[6].append(Commutator(B, Commutator(B, Commutator(B, Commutator(B, Commutator(B, A))))))

    E.append([0])
    for i, j in itertools.product(range(1, len(E[6])), range(1, len(E[1]))):
        E[7].append(Commutator(E[1][j], E[6][i]))

    for n in range(8, n_max + 1):
        E.append([0])
        print("Base n = " + str(n - 1) + " -> " + str(len(E[n - 1]) - 1) + " elements.")
        for i, j in itertools.product(range(1, len(E[n - 1])), range(1, len(E[1]))):
            print("Construint la base n = " + str(n) + ". Analitzant candidat: [E_{1, " + str(j) + "}, E_{" + str(n - 1) + ", " + str(i) + "}]")
            can = Commutator(E[1][j], E[n - 1][i])
            if not escl(can, E[n]):
                E[n].append(can)
               
    pcorxets(E)
          
    c = [0]
    for i in range(1, len(E)):
        c.append(sp.symbols(chr(96 + i) + "1:" + str(len(E[i]))))

    for n in range(5, n_max + 1):
        print("Ordre " + str(n) + ":")
        cesq = sp.S(0)
        for i in range(1, len(E[n])):
            cesq += c[n][i - 1]*E[n][i]
        cesq = cesq.doit().expand()
        llistes = []
        for i in range(1, (n//2) + 1):
            llistes.append([i, n - i])
        for ind in (llistes):
            i, k = ind
            for j, l in itertools.product(range(1, len(E[i])), range(1, len(E[k]))):
                if not (i == k and j > l):
                    cdre = Commutator(E[i][j], E[k][l]).doit().expand()
                    sol = resoldre(cesq, cdre)
                    if sol:
                        cad = "[E" + str(i) + str(j) + ", E" + str(k) + str(l) + "] = "
                        cadrel = "relE[" + str(i) + "][" + str(j) + "][" + str(k) + "][" + str(l) + "] = ["
                        for m in sol:
                            if sol[m] != 0:
                                if sol[m] < 0:
                                    cad += str(sol[m]) + "*E" + str(n) + str(c[n].index(m) + 1)
                                else:
                                    cad += "+" + str(sol[m]) + "*E" + str(n) + str(c[n].index(m) + 1)
                                nume, deno = sp.fraction(sol[m])
                                cadrel += "[" + str(nume) + ", " + str(deno) + ", " + str(n) + ", " + str(c[n].index(m) + 1) + "], "
                        cadrel = cadrel[:-2] + "]"
                        print(cad)
                        #print(cadrel)
