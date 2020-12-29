#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 09-12-2020
# alex
# rc.py

import numpy as np
import sympy as sp
import time as tm
import mpmath
import sys
sys.path.insert(0, '../')
from bchpy import *
import relations as rl

if __name__ == "__main__":
    a = sp.symbols("a0:4", real=True)
    b = sp.symbols("b0:4")
    mpmath.mp.dps = 50
    
    # print("RKN: Ordre 6 projectat, ordre 3 sense projectar. rknc36")
    # nom = "rknc36ABA13"
    # t0 = tm.time()
    # met = Metode()
    # met.setABA(a[0], b[0].expand(complex=True), a[1], b[1].expand(complex=True), a[2], b[2].expand(complex=True), a[3],
    #            b[2].conjugate().expand(complex=True), a[2], b[1].conjugate().expand(complex=True), a[1], b[0].conjugate().expand(complex=True), a[0],
    #            debug = True)
    # met.collectI()
    # met.exportpy(nom + ".py")
    # t1 = tm.time()
    # print("\tComposició ABA: 13 fluxes ->", nom, "\t(" + str(t1-t0) + " s)") 
    # A = sp.nsolve((met.w[1][1].subs({a[3]:0.1}) - 1, met.w[1][2].subs({a[3]:0.1}) - 1, met.w[2][1].subs({a[3]:0.1}), met.w[3][1].subs({a[3]:0.1}), met.w[3][2].subs({a[3]:0.1}), met.w[5][1].subs({a[3]:0.1}), met.w[5][2].subs({a[3]:0.1}), met.w[5][3].subs({a[3]:0.1}), met.w[5][4].subs({a[3]:0.1})), (a[0], a[1], a[2], re(b[0]), im(b[0]), re(b[1]), im(b[1]), re(b[2]), im(b[2])), (0.10484227087758856, 0.07919999262986678, 0.26595773649254467, 0.2439523773350367, -0.13191937900572379, 0.03013968250054609, 0.19570179373295757, 0.2259079401644172, -0.1940906708013141))
    # print(A)
    
    # nom = "rknc36BAB11"
    # t0 = tm.time()
    # met = Metode()
    # met.setBAB(b[0].expand(complex=True), a[0], b[1].expand(complex=True), a[1], b[2].expand(complex=True), a[2],
    #            b[2].conjugate().expand(complex=True), a[1], b[1].conjugate().expand(complex=True), a[0], b[0].conjugate().expand(complex=True),
    #            debug = True)   
    # met.collectI()
    # met.exportpy(nom + ".py")
    # t1 = tm.time()
    # print("\tComposició BAB: 11 fluxes ->", nom, "\t(" + str(t1-t0) + " s)") 
    
    # nom = "rknc36BAB13"
    # t0 = tm.time()
    # met = Metode()
    # met.setBAB(b[0].expand(complex=True), a[0], b[1].expand(complex=True), a[1], b[2].expand(complex=True), a[2], b[3],
    #            a[2], b[2].conjugate().expand(complex=True), a[1], b[1].conjugate().expand(complex=True), a[0], b[0].conjugate().expand(complex=True),
    #            debug = True)
    # met.collectI()
    # met.exportpy(nom + ".py")
    # t1 = tm.time()
    # print("\tComposició BAB: 13 fluxes ->", nom, "\t(" + str(t1-t0) + " s)")  

    # print("RKN: Ordre 6 projectat, ordre 5 sense projectar. rknc56")
    # nom = "rknc56BAB15"
    # t0 = tm.time()
    # met = Metode()
    # met.setBAB(b[0].expand(complex=True), a[0], b[1].expand(complex=True), a[1], b[2].expand(complex=True), a[2], b[3].expand(complex=True), a[3],
    #            b[3].conjugate().expand(complex=True), a[2], b[2].conjugate().expand(complex=True), a[1], b[1].conjugate().expand(complex=True), a[0], b[0].conjugate().expand(complex=True),
    #            debug = True)
    # met.collectI()
    # met.exportpy(nom + ".py")
    # t1 = tm.time()
    # print("\tComposició BAB: 15 fluxes ->", nom, "\t(" + str(t1-t0) + " s)")  
    # A = sp.nsolve((met.w[1][1].subs({a[3]:0.035}) - 1, met.w[1][2].subs({a[3]:0.035}) - 1, met.w[2][1].subs({a[3]:0.035}), met.w[3][1].subs({a[3]:0.035}), met.w[3][2].subs({a[3]:0.035}), met.w[4][1].subs({a[3]:0.035}), met.w[4][2].subs({a[3]:0.035}), met.w[5][1].subs({a[3]:0.035}), met.w[5][2].subs({a[3]:0.035}), met.w[5][3].subs({a[3]:0.035}), met.w[5][4].subs({a[3]:0.035})), (a[0], a[1], a[2], re(b[0]), im(b[0]), re(b[1]), im(b[1]), re(b[2]), im(b[2]), re(b[3]), im(b[3])), (0.012782593964494137, 0.20248506255393325, 0.26723234348157265, 0.004060163049474615, 0.2313745557454708, 0.07287138585245682, -0.25035588633716177, 0.28754326381740736, 0.001366669660332097, 0.1355251872806612, 0.33721900564317914))
    # print(A)
    
    # nom = "rknc56ABA15"
    # t0 = tm.time()
    # met = Metode()
    # met.setABA(a[0], b[0].expand(complex=True), a[1], b[1].expand(complex=True), a[2], b[2].expand(complex=True), a[3], b[3],
    #            a[3], b[2].conjugate().expand(complex=True), a[2], b[1].conjugate().expand(complex=True), a[1], b[0].conjugate().expand(complex=True), a[0],
    #            debug = True)
    # met.collectI()
    # # met.exportpy(nom + ".py")
    # t1 = tm.time()
    # print("\tComposició ABA: 15 fluxes ->", nom, "\t(" + str(t1-t0) + " s)")
    # A = sp.nsolve((met.w[1][1] - 1, met.w[1][2] - 1, met.w[2][1], met.w[3][1], met.w[3][2], met.w[4][1], met.w[4][2], met.w[5][1], met.w[5][2], met.w[5][3], met.w[5][4]), (a[0], a[1], a[2], a[3], re(b[0]), im(b[0]), re(b[1]), im(b[1]), re(b[2]), im(b[2]), b[3]), (0.06891288552757167, 0.17835175225535496, 0.02892351037374349, 0.22381185184332988, 0.16117063688769717, 0.014575314911348195, 0.11031233692911584, -0.2448396892914827, 0.09321147608536293, 0.2484069393674099, 0.27061110019564816))
    # print(A)
    
    # print("RKN: Ordre 4 projectat, ordre 3 sense projectar. rknc34")
    # nom = "rknc34ABA7"
    # t0 = tm.time()
    # met = Metode()
    # met.setABA(a[0], b[0].expand(complex=True), a[1], b[1],
    #            a[1], b[0].conjugate().expand(complex=True), a[0], debug = True)
    # met.collectI()
    # met.exportpy(nom + ".py")
    # met.cprint()
    # t1 = tm.time()
    # print("\tComposició ABA: 7 fluxes ->", nom, "\t(" + str(t1-t0) + " s)")  
    
    # nom = "rknc34BAB7"
    # t0 = tm.time()
    # met = Metode()
    # met.setBAB(b[0].expand(complex=True), a[0], b[1].expand(complex=True), a[1],
    #         b[1].conjugate().expand(complex=True), a[0], b[0].conjugate().expand(complex=True))
    # met.collectI()
    # met.exportpy(nom + ".py")
    # t1 = tm.time()
    # print("\tComposició BAB: 7 fluxes ->", nom, "\t(" + str(t1-t0) + " s)")
    # A = sp.nsolve((met.w[1][1].subs({a[0]:0.3725}) - 1, met.w[1][2].subs({a[0]:0.3725}) - 1, met.w[2][1].subs({a[0]:0.3725}), met.w[3][1].subs({a[0]:0.3725}), met.w[3][2].subs({a[0]:0.3725})), (a[1], re(b[0]), im(b[0]), re(b[1]), im(b[1])), (0.255, 0.143484, -0.062353, 0.356515, 0.244521))
    # print(A)
    
    # print("RKN: Ordre 4 projectat, ordre 4 sense projectar. rknc44")
    # nom = "rknc44ABA9"
    # t0 = tm.time()
    # met = Metode()
    # met.setABA(a[0], b[0].expand(complex=True), a[1], b[1].expand(complex=True), a[2],
    #            b[1].conjugate().expand(complex=True), a[1], b[0].conjugate().expand(complex=True), a[0])
    # met.collectI()
    # met.exportpy(nom + ".py")
    # t1 = tm.time()
    # print("\tComposició ABA: 9 fluxes ->", nom, "\t(" + str(t1-t0) + " s)") 

    # nom = "rknc44ABA11"
    # t0 = tm.time()
    # met = Metode()
    # met.setABA(a[0], b[0].expand(complex=True), a[1], b[1].expand(complex=True), a[2], b[2],
    #            a[2], b[1].conjugate().expand(complex=True), a[1], b[0].conjugate().expand(complex=True), a[0])
    # met.collectI()
    # met.exportpy(nom + ".py")
    # t1 = tm.time()
    # print("\tComposició ABA: 11 fluxes ->", nom, "\t(" + str(t1-t0) + " s)")     

    # nom = "rknc44BAB9"
    # t0 = tm.time()
    # met = Metode()
    # met.setBAB(b[0].expand(complex=True), a[0], b[1].expand(complex=True), a[1], b[2],
    #            a[1], b[1].conjugate().expand(complex=True), a[0], b[0].conjugate().expand(complex=True))
    # met.collectI()
    # met.exportpy(nom + ".py")
    # t1 = tm.time()
    # print("\tComposició BAB: 9 fluxes ->", nom, "\t(" + str(t1-t0) + " s)")  
