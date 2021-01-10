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
    impfits = True
    n_max = 10
    
    A = Operator("A")
    B = Operator("B")
    
    E = []
    E.append([0])
    Ep = [[0]]
    
    E.append([0])
    E[1].append(A)
    E[1].append(B)
    Ep.append([0, "A", "B"])

    
    E.append([0])
    E[2].append(Commutator(A, B))
    Ep.append([0, "AB"])
    
    E.append([0])
    E[3].append(Commutator(A, E[2][1]))
    E[3].append(Commutator(B, E[2][1]))
    Ep.append([0, "AAB", "BAB"])

    
    E.append([0])
    E[4].append(Commutator(A, Commutator(A, Commutator(A, B))))
    E[4].append(Commutator(B, Commutator(A, Commutator(A, B))))
    E[4].append(Commutator(B, Commutator(B, Commutator(B, A))))
    Ep.append([0, "AAAB", "BAAB", "BBBA"])

    
    E.append([0])
    E[5].append(Commutator(A, Commutator(A, Commutator(A, Commutator(A, B)))))
    E[5].append(Commutator(B, Commutator(A, Commutator(A, Commutator(A, B)))))
    E[5].append(Commutator(A, Commutator(A, Commutator(B, Commutator(B, A)))))
    E[5].append(Commutator(B, Commutator(B, Commutator(A, Commutator(A, B)))))
    E[5].append(Commutator(A, Commutator(B, Commutator(B, Commutator(B, A)))))
    E[5].append(Commutator(B, Commutator(B, Commutator(B, Commutator(B, A)))))
    Ep.append([0, "AAAAB", "BAAAB", "AABBA", "BBAAB", "ABBBA", "BBBBA"])

    
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
    Ep.append([0, "AAAAAB", "BAAAAB", "ABAAAB", "ABBAAB", "BBAAAB", "AABBBA", "BABBBA", "ABBBBA", "BBBBBA"])

    
    E.append([0])
    Ep.append([0])
    for i, j in itertools.product(range(1, len(E[6])), range(1, len(E[1]))):
        E[7].append(Commutator(E[1][j], E[6][i]))
        Ep[7].append(("A" if j == 1 else "B") + Ep[6][i])


    aux1 = []
    aux2 = []
    for i in range(0, len(E) - 1):
        aux1.append([[]]*len(E[i]))
        aux2.append([[]]*len(E[i]))
    solvec = [0, aux1, aux2]
    
    for n in range(8, n_max + 1):
        E.append([0])
        Ep.append([0])
        print("Base n = " + str(n - 1) + " -> " + str(len(E[n - 1]) - 1) + " elements.")
        solvec[1].append([[]]*len(E[n - 1]))
        solvec[2].append([[]]*len(E[n - 1]))
        l = 1
        for i, j in itertools.product(range(1, len(E[n - 1])), range(1, len(E[1]))):
            print("Construint la base n = " + str(n) + ". Analitzant candidat: [E_{1, " + str(j) + "}, E_{" + str(n - 1) + ", " + str(i) + "}]")
            can = Commutator(E[1][j], E[n - 1][i])
            sol = escl(can, E[n], A)
            if not escl(can, E[n], A):
                E[n].append(can)
                Ep[n].append(("A" if j == 1 else "B") + Ep[n - 1][i])
                solvec[j][n - 1][i] = {sp.Symbol("s" + str(l)):1}
                l += 1
            else:
                solvec[j][n - 1][i] = sol
    pcorxets(E)

    if impfits:
        fitC = open("cniats.txt", "w")
        fitL = open("tex.txt", "w") 
        fitR = open("relArray.txt", "w")
        for i in range(1, len(Ep)):
            fitC.write(str(i) + "\n")
            for j in range(1, len(Ep[i])):
                fitC.write(Ep[i][j] + "\n")
        fitC.close()
        for i in range(1, len(Ep)):
            fitL.write(str(i) + "\n")
            for j in range(1, len(Ep[i])):
                cad = "E_{" + str(i) + "," + str(j) + "}=" + "\left["
                for k in range(0, len(Ep[i][j])):
                    cad += Ep[i][j][k] + ","
                cad = cad[:-1] + "\\right]"
                fitL.write(cad + "\n")       

    c = [0]
    for i in range(1, len(E)):
        c.append(sp.symbols("s1:" + str(len(E[i]))))

    for n in range(5, n_max + 1):
        print("Ordre " + str(n) + ":")
        llistes = []
        for i in range(1, (n//2) + 1):
            llistes.append([i, n - i])
        for ind in (llistes):
            i, k = ind
            for j, l in itertools.product(range(1, len(E[i])), range(1, len(E[k]))):
                if not (i == k and j > l):
                    if i == 1 and solvec[j][k][l]:
                        sol = solvec[j][k][l]
                    else:
                        cdre = Commutator(E[i][j], E[k][l])
                        cesq = sp.S(0)
                        for m in range(1, len(E[n])):
                            if cdre.count(A) == E[n][m].count(A):
                                cesq += c[n][m - 1]*E[n][m]
                        cesq = cesq.doit().expand()
                        cdre = cdre.doit().expand()
                        sol = resoldre(cesq, cdre)
                    if sol:
                        cad = "[E" + str(i) + str(j) + ", E" + str(k) + str(l) + "] = "
                        cadrel = "relE[" + str(i) + "][" + str(j) + "][" + str(k) + "][" + str(l) + "] = ["
                        cadtex = "\left[E_{" + str(i) + "," + str(j) + "},E_{" + str(k) + "," + str(l) + "} \\right] = "
                        for m in sol:
                            if sol[m] != 0:
                                if sol[m] < 0:
                                    cad += str(sol[m]) + "*E" + str(n) + str(c[n].index(m) + 1)
                                else:
                                    cad += "+" + str(sol[m]) + "*E" + str(n) + str(c[n].index(m) + 1)
                                nume, deno = sp.fraction(sol[m])
                                cadrel += "[" + str(nume) + ", " + str(deno) + ", " + str(n) + ", " + str(c[n].index(m) + 1) + "], "
                                factex = ""
                                if abs(nume) != 1 and deno != 1:
                                    factex = str(nume)
                                    if deno != 1:
                                        factex = "\\frac{" + factex + "}{" + str(deno) + "}"
                                cadtex += "+" + factex + "E_{" + str(n) + "," + str(c[n].index(m) + 1) + "}"
                        cadrel = cadrel[:-2] + "]"
                        print(cad)
                        if impfits:
                            fitR.write(cadrel + "\n")
                            fitL.write(cadtex + "\n")
    if impfits:
        fitR.close()
        fitL.close()
