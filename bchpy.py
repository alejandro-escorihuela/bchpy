#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 01-12-2020
# alex
# bchpy.py

import numpy as np
import sympy as sp

class Eel():
    def __init__(self, i = 1, j = 1):
        self.i = i
        self.j = j
        self.simbol = sp.symbols("E" + str(i) + str(j))
        
    def imprimir(self):
        print(self.simbol)
        
if __name__ == "__main__":
    E11=Eel(3,1)
    E11.imprimir()

