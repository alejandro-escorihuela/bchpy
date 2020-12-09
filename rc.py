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
    x, h = sp.symbols("x h")
    a = sp.symbols("a0:2", real=True)
    b = sp.symbols("b0:1")

    t0 = tm.time()
    w = mat_esqABA(a[0], b[0].expand(complex=True), a[1], b[0].conjugate().expand(complex=True), a[0])
    for i in range(1, len(w)):
        for j in range(1, len(w[i])):
            print("w(" + str(i) + "," + str(j) + ") =", w[i][j])
    t1 = tm.time()
    print(t1-t0, "s")    
