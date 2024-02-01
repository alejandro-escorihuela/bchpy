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

    ## Mètode d'escissió RKN ABA d'ordre 8
    cofs = [0.052092434384033901862, 0.14585030481264474322, 0.22528749326770217132, 0.2551565441392939504,
            0.41627618961225709704,  0.018133468820831725316, -0.38456727021395037402, -0.17904011029926455989,
            0.09972717834705148443, -0.11847080143330224189, -0.108833834399100215506, 0.1864616892738210907,
            0.22201073664899168003, 0.45904158176713683037, 0.52387952203673426865, -0.0036608362703183590543,
            -0.5458724496837200138, -0.52694368162168635835, -0.5458724496837200138, -0.0036608362703183590543,
            0.52387952203673426865, 0.45904158176713683037, 0.22201073664899168003, 0.1864616892738210907,
            -0.108833834399100215506, -0.11847080143330224189, 0.09972717834705148443, -0.17904011029926455989,
            -0.38456727021395037402, 0.018133468820831725316, 0.41627618961225709704, 0.2551565441392939504,
            0.22528749326770217132, 0.14585030481264474322, 0.052092434384033901862]
    met = Metode(depth = 9, basis_type = "E", numeric = True, rkn = True)
    t0 = tm.time()
    met.setABA(*cofs)
    print("Temps = %f" % (tm.time() - t0))
    correcte = abs(met.w[1][1] + met.w[1][2] - 2.0) < tol
    for i in range(2, 7):
        for j in range(1, len(met.w[i])):
            correcte = correcte and abs(met.w[i][j]) < tol
    if not correcte:
        printe("Resultat erroni en Mètode RKN ABA simètric d'ordre 8")
        met.cprint()
        
    ## Mètode d'escissió ABA d'ordre 6
    cofs = [0.0502627644003922, 0.148816447901042, 0.413514300428344, -0.132385865767784, 0.0450798897943977, 0.067307604692185, -0.188054853819569, 0.432666402578175, 0.54196067845078, -0.016404589403617997, -0.7255255585086897, -0.016404589403617997, 0.54196067845078, 0.432666402578175, -0.188054853819569, 0.067307604692185, 0.0450798897943977, -0.132385865767784, 0.413514300428344, 0.148816447901042, 0.0502627644003922]
    met = Metode(depth = 7, basis_type = "E", numeric = True)
    t0 = tm.time()
    met.setABA(*cofs)
    print("Temps = %f" % (tm.time() - t0))
    correcte = abs(met.w[1][1] + met.w[1][2] - 2.0) < tol
    for i in range(2, 7):
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
    for i in range(2, 7):
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
    for i in range(2, 7):
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
    for i in range(2, 7):
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
    for i in range(2, 9):
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
    for i in range(2, 7):
        for j in range(1, len(met.w[i])):
            correcte = correcte and abs(met.w[i][j]) < tol
    if not correcte:
        printe("Resultat erroni en Mètode SS(S4) simètric d'ordre 6")
        met.cprint()

    ## Mètode SC(S2) simètric-conjugat d'ordre 5
    cofs = [(0.1752684090720741140583563+0.05761474413053870201304364j), (0.1848736801929841604288898-0.1941219227572495885067758j), 0.2797158214698834510255077, (0.1848736801929841604288898+0.1941219227572495885067758j), (0.1752684090720741140583563-0.05761474413053870201304364j)] 
    met = Metode(depth = 6, basis_type = "M", numeric = True)
    met.setSS2(*cofs)
    correcte = abs(met.w[1][1].real - 1.0) < tol and abs(met.w[1][1].imag) < tol
    for i in range(2, 6):
        for j in range(1, len(met.w[i])):
            correcte = correcte and abs(met.w[i][j].real) < tol and abs(met.w[i][j].imag) < tol
    if not correcte:
        printe("Resultat erroni en Mètode SC(S2) simètric-conjugat d'ordre 5")
        met.cprint()        

    ## Mètode SS(S2) simètric d'ordre 6 amb coeficients complexos
    cofs =[(0.11690003755466129+0.04342825461606034j), (0.12955910128208825-0.1239896121880926j), (0.18653249281213383+0.003107430710072675j), (0.13401673670223335+0.15490785372391916j), (0.18653249281213383+0.003107430710072675j), (0.12955910128208825-0.1239896121880926j), (0.11690003755466129+0.04342825461606034j)]
    met = Metode(depth = 7, basis_type = "M", numeric = True)
    met.setSS2(*cofs)
    correcte = abs(met.w[1][1].real - 1.0) < tol and abs(met.w[1][1].imag) < tol
    for i in range(2, 7):
        for j in range(1, len(met.w[i])):
            correcte = correcte and abs(met.w[i][j].real) < tol and abs(met.w[i][j].imag) < tol
    if not correcte:
        printe("Resultat erroni en Mètode SS(S2) simètric d'ordre 6 amb coeficients complexos")
        met.cprint()      
