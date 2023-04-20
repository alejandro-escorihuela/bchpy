#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 17-04-2023
# alex
# compS4.py
# Generar la base per a composicions de mètodes S4.

import sympy as sp
from sympy.physics.quantum import Operator, Commutator
import itertools
import sys
sys.path.insert(0, '../')
from qolib import *

if __name__ == "__main__":
    impfits = True
    n_max = 15

    S = []
    S.append([0])
    for i in range(0, n_max + 1):
        S.append(Operator("S_{%d}" % (i + 1)))

    Z = []
    Z.append([0])
    
    Z.append([0])
    Z[1].append(S[1])
    
    Z.append([0])
    
    Z.append([0])
    #Z[3].append(S[3])

    for n in range(4, n_max + 1):
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
        fitC = open("compS4_base.txt", "w")
        for i in range(1, len(Z)):
            for j in range(1, len(Z[i])):
                fitC.write(str(Z[i][j]) + " ")
            fitC.write("\n")
        fitC.close()     
