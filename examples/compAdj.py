#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 24-04-2021
# alex
# compAdj.py
# Generar la base per a composicions de primer ordre.

import sympy as sp
from sympy.physics.quantum import Operator, Commutator
import itertools
import sys
sys.path.insert(0, '../')
from qolib import *

if __name__ == "__main__":
    impfits = True
    n_max = 9

    Y = []
    Y.append([0])
    for i in range(n_max + 1):
        Y.append(Operator("Y_{%d}" % (i + 1)))

    Z = []
    Z.append([0])
    
    Z.append([0])
    Z[1].append(Y[1])
    
    Z.append([0])
    Z[2].append(Y[2])
    
    Z.append([0])
    Z[3].append(Y[3])
    Z[3].append(Commutator(Z[1][1], Z[2][1]))

    for n in range(4, n_max + 1):
        Z.append([0])
        Z[n].append(Y[n])
        llistes = []
        for i in range(1, (n//2) + 1):
            llistes.append([i, n - i])
        for ind in (llistes):
            i, k = ind
            for j, l in itertools.product(range(1, len(Z[i])), range(1, len(Z[k]))):
                if not (i == k and j >= l):
                    can = Commutator(Z[i][j], Z[k][l])
                    if not escl(can, Z[n]):
                        Z[n].append(can)
    pcorxets(Z, label = "Z")

    if impfits:
        fitC = open("compAdj_base.txt", "w")
        for i in range(1, len(Z)):
            for j in range(1, len(Z[i])):
                fitC.write(str(Z[i][j]) + " ")
            fitC.write("\n")
        fitC.close()     
