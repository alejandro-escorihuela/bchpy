#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 20-04-2023
# alex
# compS10Id.py
# Identitats (relacions de commutació) dels corxets de la base per a composicions de S10

import sympy as sp
from sympy.physics.quantum import Operator, Commutator
import itertools
import sys
sys.path.insert(0, '../')
from qolib import *

if __name__ == "__main__":
    impfits = True
    n_max = 27

    S = []
    S.append([0])
    for i in range(n_max + 1):
        S.append(Operator("S_{%d}" % (i + 1)))
        
    Z = llegir_baseS2("compS10_base.txt", S)
    
    if (len(Z) - 1 < n_max):
        print("La dimensió de la base llegida és menor que l'ordre a calcular. (" + str(len(Z) - 1) + "<" + str(n_max) + ")")
    
    if impfits:
        fitR = open("relArray.txt", "w")     

    c = [0]
    for i in range(1, len(Z)):
        c.append(sp.symbols("s1:" + str(len(Z[i]))))

    for n in range(1, n_max + 1):
        print("Ordre " + str(n) + ":")
        llistes = []
        for i in range(1, (n//2) + 1):
            llistes.append([i, n - i])
        for ind in (llistes):
            i, k = ind
            for j, l in itertools.product(range(1, len(Z[i])), range(1, len(Z[k]))):
                if not (i == k and j > l):
                    cdre = Commutator(Z[i][j], Z[k][l])
                    cesq = sp.S(0)
                    for m in range(1, len(Z[n])):
                        cesq += c[n][m - 1]*Z[n][m]
                    cesq = cesq.doit().expand()
                    cdre = cdre.doit().expand()
                    sol = resoldre(cesq, cdre)
                    if sol:
                        cad = "[M10" + str(i) + str(j) + ", M10" + str(k) + str(l) + "] = "
                        cadrel = "relM10[" + str(i) + "][" + str(j) + "][" + str(k) + "][" + str(l) + "] = ["
                        # cadtex = "\left[M10_{" + str(i) + "," + str(j) + "},M10_{" + str(k) + "," + str(l) + "} \\right] = "
                        for m in sol:
                            if sol[m] != 0:
                                if sol[m] < 0:
                                    cad += str(sol[m]) + "*M10" + str(n) + str(c[n].index(m) + 1)
                                else:
                                    cad += "+" + str(sol[m]) + "*M10" + str(n) + str(c[n].index(m) + 1)
                                nume, deno = sp.fraction(sol[m])
                                cadrel += "[" + str(nume) + ", " + str(deno) + ", " + str(n) + ", " + str(c[n].index(m) + 1) + "], "
                                # factex = ""
                                # if abs(nume) != 1 and deno != 1:
                                #     factex = str(nume)
                                #     if deno != 1:
                                #         factex = "\\frac{" + factex + "}{" + str(deno) + "}"
                                # cadtex += "+" + factex + "M10_{" + str(n) + "," + str(c[n].index(m) + 1) + "}"
                        cadrel = cadrel[:-2] + "]"
                        print(cad)
                        if impfits:
                            fitR.write(cadrel + "\n")
    if impfits:
        fitR.close()
