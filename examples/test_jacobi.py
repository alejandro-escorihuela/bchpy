#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 15-01-2021
# alex
# test_jacobi.py

import numpy as np
import sympy as sp
import time as tm
import sys
sys.path.insert(0, '../')
from bchpy import *
import relations as rl

if __name__ == "__main__":
    ord_bch = len(rl.tamN) + 1
    for i in range(1, ord_bch):
        for j in range(1, rl.tamN[i - 1] + 1):
            for k in range(1, ord_bch):
                for l in range(1, rl.tamN[k - 1] + 1):
                    X = Eel(i, j, "N")
                    Y = Eel(k, l, "N")
                    Z1 = (Claudator(Y, Claudator(X, Claudator(X, Y)))).expand()
                    Z2 = (Claudator(X, Claudator(Y, Claudator(X, Y)))).expand()
                    Z = (Z1-Z2).expand()
                    if Z != 0:
                        print(X, Y, "->",  (Z1-Z2).expand())
