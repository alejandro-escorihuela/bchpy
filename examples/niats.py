#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 28-12-2020
# alex
# niats.py
# Generar la base de corxets niats o niuats

import sympy as sp
from sympy.physics.quantum import Operator, Commutator
import itertools
import sys
sys.path.insert(0, '../')
from qolib import *

if __name__ == "__main__":
    impfits = True
    n_max = 11
    
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
  
    for n in range(8, n_max + 1):
        E.append([0])
        Ep.append([0])
        print("Base n = " + str(n - 1) + " -> " + str(len(E[n - 1]) - 1) + " elements.")
        for i, j in itertools.product(range(1, len(E[n - 1])), range(1, len(E[1]))):
            print("Construint la base n = " + str(n) + ". Analitzant candidat: [E_{1, " + str(j) + "}, E_{" + str(n - 1) + ", " + str(i) + "}]")
            can = Commutator(E[1][j], E[n - 1][i])
            sol = esclCA(can, E[n], A)
            if not esclCA(can, E[n], A):
                E[n].append(can)
                Ep[n].append(("A" if j == 1 else "B") + Ep[n - 1][i])
    pcorxets(E)

    if impfits:
        fitC = open("niats_base.txt", "w")
        for i in range(1, len(Ep)):
            for j in range(1, len(Ep[i])):
                fitC.write("[" + Ep[i][j] + "]" + " ")
            fitC.write("\n")
        fitC.close()     


