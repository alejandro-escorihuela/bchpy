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
    print("rc")
    a = sp.symbols("a0:4", real=True)
    b = sp.symbols("b0:3")

    t0 = tm.time()
    #w = mat_esqABA(a[0], b[0].expand(complex=True), a[1], b[1].expand(complex=True), a[2], b[2].expand(complex=True), a[3], b[2].conjugate().expand(complex=True), a[2], b[1].conjugate().expand(complex=True), a[1], b[0].conjugate().expand(complex=True), a[0], debug = True)
    w = mat_esqABA(a[0], b[0].expand(complex=True), a[1], b[0].conjugate().expand(complex=True), a[0], debug = True)
    w = mat_collectI(w)
    print_mat(w)
    mat_exportC(w, "prova.py")
    t1 = tm.time()
    print(t1-t0, "s")    
