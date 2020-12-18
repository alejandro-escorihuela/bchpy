#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 09-12-2020
# alex
# rc.py

import numpy as np
import sympy as sp
import relations as rl
import time as tm
from bchpy import *

if __name__ == "__main__":
    a = sp.symbols("a0:4", real=True)
    b = sp.symbols("b0:4")

    print("RKN: Ordre 6 projectat, ordre 3 sense projectar. rknc36")
    nom = "rknc36ABA13"
    t0 = tm.time()
    met = Metode()
    met.setABA(a[0], b[0].expand(complex=True), a[1], b[1].expand(complex=True), a[2], b[2].expand(complex=True), a[3],
               b[2].conjugate().expand(complex=True), a[2], b[1].conjugate().expand(complex=True), a[1], b[0].conjugate().expand(complex=True), a[0], debug = True)
    met.collectI()
    met.exportpy(nom + ".py")
    t1 = tm.time()
    print("\tComposició ABA: 13 fluxes ->", nom, "\t(" + str(t1-t0) + " s)") 

    nom = "rknc36BAB13"
    t0 = tm.time()
    met = Metode()
    met.setBAB(b[0].expand(complex=True), a[0], b[1].expand(complex=True), a[1], b[2].expand(complex=True), a[2],
               b[2].conjugate().expand(complex=True), a[1], b[1].conjugate().expand(complex=True), a[0], b[0].conjugate().expand(complex=True), debug = True)   
    met.collectI()
    met.exportpy(nom + ".py")
    t1 = tm.time()
    print("\tComposició BAB: 13 fluxes ->", nom, "\t(" + str(t1-t0) + " s)") 
    
    nom = "rknc36BAB11"
    t0 = tm.time()
    met = Metode()
    met.setBAB(b[0].expand(complex=True), a[0], b[1].expand(complex=True), a[1], b[2].expand(complex=True), a[2], b[3],
               a[2], b[2].conjugate().expand(complex=True), a[1], b[1].conjugate().expand(complex=True), a[0], b[0].conjugate().expand(complex=True), debug = True)
    met.collectI()
    met.exportpy(nom + ".py")
    t1 = tm.time()
    print("\tComposició BAB: 11 fluxes ->", nom, "\t(" + str(t1-t0) + " s)")  

    # print("RKN: Ordre 6 projectat, ordre 5 sense projectar. rknc36")
    # nom = "rknc56BAB15"
    # t0 = tm.time()
    # met = Metode()
    # met.setBAB(b[0].expand(complex=True), a[0], b[1].expand(complex=True), a[1], b[2].expand(complex=True), a[2], b[3].expand(complex=True), a[3],
    #         b[3].conjugate().expand(complex=True), a[2], b[2].conjugate().expand(complex=True), a[1], b[1].conjugate().expand(complex=True), a[0], b[0].conjugate().expand(complex=True))
    # met.collectI()
    # met.exportpy(nom + ".py")
    # t1 = tm.time()
    # print("\tComposició BAB: 15 fluxes ->", nom, "\t(" + str(t1-t0) + " s)")  

    # print("RKN: Ordre 4 projectat, ordre 3 sense projectar. rknc34")
    # nom = "rknc34BAB7"
    # t0 = tm.time()
    # met = Metode()
    # met.setBAB(b[0].expand(complex=True), a[0], b[1].expand(complex=True), a[1],
    #         b[1].conjugate().expand(complex=True), a[0], b[0].conjugate().expand(complex=True))
    # met.collectI()
    # met.exportpy(nom + ".py")
    # t1 = tm.time()
    # print("\tComposició BAB: 7 fluxes ->", nom, "\t(" + str(t1-t0) + " s)")  
    

    # print("RKN: Ordre 4 projectat, ordre 4 sense projectar. rknc44")
    # nom = "rknc44ABA9"
    # t0 = tm.time()
    # met = Metode()
    # met.setABA(a[0], b[0].expand(complex=True), a[1], b[1].expand(complex=True), a[2], b[2].expand(complex=True), a[3],
    #            b[2].conjugate().expand(complex=True), a[2], b[1].conjugate().expand(complex=True), a[1], b[0].conjugate().expand(complex=True), a[0])
    # met.collectI()
    # met.exportpy(nom + ".py")
    # t1 = tm.time()
    # print("\tComposició ABA: 9 fluxes ->", nom, "\t(" + str(t1-t0) + " s)") 
