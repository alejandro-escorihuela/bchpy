#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 28-12-2020
# alex
# relE7.py

import sympy as sp
from sympy.physics.quantum import Operator, Commutator
import itertools
# import sys
# sys.path.insert(0, '../')
# from bchpy import *
        
def resoldre(exprEsq, exprDre, debug = False):
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
        if equ not in eqs and -equ not in eqs:
            eqs.append(equ) 
    if debug == True:
        for i in eqs:
            print(i)
    return sp.solve(eqs)

def escl(can, En):
    if len(En) == 1:
        return False
    sim = sp.symbols("s1:" + str(len(En)))
    cesq = sp.S(0)
    for i in range(1, len(En)):
        cesq += sim[i - 1]*En[i]
    cesq = cesq.doit().expand()
    cane = can.doit().expand()
    if resoldre(cesq, cane):
        return True
    return False

def pcorxets(E):
    for i in range(1, len(E)):
        for j in range(1, len(E[i])):
            print("E_{" + str(i) + ", " + str(j) + "} =", E[i][j])
            
if __name__ == "__main__":
    n_max = 8
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
        for i, j in itertools.product(range(1, len(E[n - 1])), range(1, len(E[1]))):
            #print("Construint la base n = " + str(n) + ". Analitzant candidat: [E_{1, " + str(j) + "}, E_{" + str(n - 1) + ", " + str(i) + "}]")
            can = Commutator(E[1][j], E[n - 1][i])
            if not escl(can, E[n]):
                E[n].append(can)


    ## Test resoldre
    c = [0]
    for i in range(1, len(E)):
        c.append(sp.symbols(chr(96 + i) + "1:" + str(len(E[i]))))
    cesq = sp.S(0)
    for i in range(1, len(E[8])):
        cesq += c[8][i - 1]*E[8][i]
    cesq = cesq.doit().expand()
    cdre = Commutator(A, E[7][11]).doit().expand()
    sol = resoldre(cesq, cdre, debug = True)
    print(sol)
    exit(-1)
    ## Fi test
                
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
                                #print(sol[m])
                                # if sol[m] < 0:
                                #     cad += str(sol[m]) + "*E" + str(n) + str(c[n].index(m) + 1)
                                # else:
                                #     cad += "+" + str(sol[m]) + "*E" + str(n) + str(c[n].index(m) + 1)
                                nume, deno = sp.fraction(sol[m])
                                cadrel += "[" + str(nume) + ", " + str(deno) + ", " + str(n) + ", " + str(c[n].index(m) + 1) + "], "
                        cadrel = cadrel[:-2] + "]"
                        #print(cad)
                        print(cadrel)
