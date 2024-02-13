#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 28-12-2020
# alex
# niatsId.py
# Identitats (relacions de commutació) dels corxets niats o niuats

import sympy as sp
from sympy.physics.quantum import Operator, Commutator
import itertools
import sys
sys.path.insert(0, '../')
from qolib import *

if __name__ == "__main__":
    impfits = True
    n_max = 9
    
    A = Operator("A")
    B = Operator("B")
    
    E = llegir_base("../txt/niats_base.txt", A, B)
    if (len(E) - 1 < n_max):
        print("La dimensió de la base llegida és menor que l'ordre a calcular. (" + str(len(E) - 1) + "<" + str(n_max) + ")")
    
    if impfits:
        fitR = open("relArray.txt", "w")     

    c = [0]
    for i in range(1, len(E)):
        c.append(sp.symbols("s1:" + str(len(E[i]))))
    
    for n in range(1, n_max + 1):
        print("Ordre " + str(n) + ":")
        llistes = []
        for i in range(1, (n//2) + 1):
            llistes.append([i, n - i])
        for ind in (llistes):
            i, k = ind
            for j, l in itertools.product(range(1, len(E[i])), range(1, len(E[k]))):
                if not (i == k and j > l):
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
                        cadtex = "\left[E_{" + str(i) + "," + str(j) + "},E_{" + str(k) + "," + str(l) + "} \\right] &= "
                        for m in sol:
                            if sol[m] != 0:
                                if sol[m] < 0:
                                    cad += str(sol[m]) + "*E" + str(n) + str(c[n].index(m) + 1)
                                else:
                                    cad += "+" + str(sol[m]) + "*E" + str(n) + str(c[n].index(m) + 1)
                                nume, deno = sp.fraction(sol[m])
                                cadrel += "[" + str(nume) + ", " + str(deno) + ", " + str(n) + ", " + str(c[n].index(m) + 1) + "], "
                                # açò està mal. no és correcte.
                                if deno != 1:
                                    factex = "\\frac{" + str(abs(nume)) + "}{" + str(deno) + "}"
                                elif abs(nume) != 1:
                                    factex = str(abs(nume))
                                else:
                                    factex = ""
                                signe = "+"
                                if nume < 0:
                                    signe = "-"
                                cadtex += signe + factex + "E_{" + str(n) + "," + str(c[n].index(m) + 1) + "}"
                        cadrel = cadrel[:-2] + "]"
                        print(cadtex)
                        if impfits:
                            fitR.write(cadrel + "\n")
    if impfits:
        fitR.close()
