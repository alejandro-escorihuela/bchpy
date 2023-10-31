#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 24-04-2023
# alex
# geneqSS.py

import numpy as np
import sympy as sp
import time as tm
import sys
sys.path.insert(0, '../')
from bchpy import *
import relations as rl

if __name__ == "__main__":
    tol = 1e-13
    
    ## Mètode d'escissió ABA d'ordre 6
    cofs = [0.0502627644003922, 0.148816447901042, 0.413514300428344, -0.132385865767784, 0.0450798897943977, 0.067307604692185, -0.188054853819569, 0.432666402578175, 0.54196067845078, -0.016404589403617997, -0.7255255585086897, -0.016404589403617997, 0.54196067845078, 0.432666402578175, -0.188054853819569, 0.067307604692185, 0.0450798897943977, -0.132385865767784, 0.413514300428344, 0.148816447901042, 0.0502627644003922]
    met = Metode(depth = 7, basis_type = "E", numeric = True)
    met.setABA(*cofs)
    correcte = abs(met.w[1][1] + met.w[1][2] - 2.0) < tol
    for i in range(2, 6):
        for j in range(1, len(met.w[i])):
            correcte = correcte and abs(met.w[i][j]) < tol
    if not correcte:
        printe("Resultat erroni en Mètode ABA simètric d'ordre 6")
        met.cprint()

    ## Mètode d'escissió BAB d'ordre 6
    cofs = [0.0502627644003922, 0.148816447901042, 0.413514300428344, -0.132385865767784, 0.0450798897943977, 0.067307604692185, -0.188054853819569, 0.432666402578175, 0.54196067845078, -0.016404589403617997, -0.7255255585086897, -0.016404589403617997, 0.54196067845078, 0.432666402578175, -0.188054853819569, 0.067307604692185, 0.0450798897943977, -0.132385865767784, 0.413514300428344, 0.148816447901042, 0.0502627644003922]
    met = Metode(depth = 7, basis_type = "E", numeric = True)
    met.setBAB(*cofs)
    correcte = abs(met.w[1][1] + met.w[1][2] - 2.0) < tol
    for i in range(2, 6):
        for j in range(1, len(met.w[i])):
            correcte = correcte and abs(met.w[i][j]) < tol
    if not correcte:
        printe("Resultat erroni en Mètode BAB simètric d'ordre 6")
        met.cprint()

    ## Mètode de composició XX d'ordre 6
    cofs = [0.0502627644003922, 0.0985536835006498, 0.31496061692769417, -0.44734648269547816, 0.49242637248987586, -0.42511876779769087, 0.23706391397812188, 0.19560248860005314, 0.34635818985072686, -0.36276277925434486, -0.36276277925434486, 0.34635818985072686, 0.19560248860005314, 0.23706391397812188, -0.42511876779769087, 0.49242637248987586, -0.44734648269547816, 0.31496061692769417, 0.0985536835006498, 0.0502627644003922]
    met = Metode(depth = 7, basis_type = "C", numeric = True)
    met.setXX(*cofs)
    correcte = abs(met.w[1][1] - 1.0) < tol
    for i in range(2, 6):
        for j in range(1, len(met.w[i])):
            correcte = correcte and abs(met.w[i][j]) < tol
    if not correcte:
        printe("Resultat erroni en Mètode XX simètric d'ordre 6")
        met.cprint()
        
    ## Mètode SS(S2) simètric d'ordre 6
    cofs = [0.78451361047755726381949763, 0.23557321335935813368479318, -1.17767998417887100694641568, 1.31518632068391121888424973]
    met = Metode(depth = 7, basis_type = "M", numeric = True)
    met.setSS2(cofs[0], cofs[1], cofs[2], cofs[3], cofs[2], cofs[1], cofs[0])  
    correcte = abs(met.w[1][1] - 1.0) < tol 
    for i in range(2, 6):
        for j in range(1, len(met.w[i])):
            correcte = correcte and abs(met.w[i][j]) < tol
    if not correcte:
        printe("Resultat erroni en Mètode SS(S2) simètric d'ordre 6")
        met.cprint()

    ## Mètode SS(S2) simètric d'ordre 8
    cofs = [0.74167036435061295344822780, -0.40910082580003159399730010, 0.19075471029623837995387626, -0.57386247111608226665638773,
            0.29906418130365592384446354, 0.33462491824529818378495798, 0.31529309239676659663205666, -0.79688793935291635291635401978884]
    cofs = cofs + cofs[-2::-1]
    met = Metode(depth = 8, basis_type = "M", numeric = True)
    met.setSS2(*cofs)
    correcte = abs(met.w[1][1] - 1.0) < tol
    for i in range(2, 8):
        for j in range(1, len(met.w[i])):
            correcte = correcte and abs(met.w[i][j]) < tol
    if not correcte:
        printe("Resultat erroni en Mètode SS(S2) simètric d'ordre 8")
        met.cprint()

    ## Mètode SS(S4) simètric d'ordre 6
    alpha = 1/(2 - 2**(1/5))
    beta = 1 -2*alpha
    cofs = [alpha, beta, alpha]
    met = Metode(depth = 7, basis_type = "M4", numeric = True)
    met.setSS4(*cofs)     
    correcte = abs(met.w[1][1] - 1.0) < tol
    for i in range(2, 6):
        for j in range(1, len(met.w[i])):
            correcte = correcte and abs(met.w[i][j]) < tol
    if not correcte:
        printe("Resultat erroni en Mètode SS(S4) simètric d'ordre 6")
        met.cprint()
