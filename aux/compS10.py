#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 20-04-2023
# alex
# compS10.py
# Generar la base per a composicions de mètodes S10.

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
    for i in range(0, n_max + 1):
        S.append(Operator("S_{%d}" % (i + 1)))

    Z = []
    Z.append([0])

    # Z1
    Z.append([0])
    Z[1].append(S[1])
    # Z2, Z3, Z4, Z5, Z6, Z7, Z8, Z9, Z10
    Z.append([0])
    Z.append([0])
    Z.append([0])
    Z.append([0])
    Z.append([0])
    Z.append([0])
    Z.append([0])
    Z.append([0])
    Z.append([0])
    
    for n in range(10, n_max + 1):
        Z.append([0])
        if n % 2 != 0:
            Z[n].append(S[n])
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
        fitC = open("compS10_base.txt", "w")
        for i in range(1, len(Z)):
            for j in range(1, len(Z[i])):
                fitC.write(str(Z[i][j]) + " ")
            fitC.write("\n")
        fitC.close()     
