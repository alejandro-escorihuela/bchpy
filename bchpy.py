#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 01-12-2020
# alex
# bchpy.py

import numpy as np
import sympy as sp
import relations as rl
from sympy.printing.pretty.stringpict import prettyForm
from sympy.printing.pretty.pretty_symbology import pretty_symbol

class Eel(sp.Expr):
    is_commutative = False
    def __init__(self, i = 1, j = 1):
        self.i = i
        self.j = j
        
    def _sympystr(self, printer, *args):
        return "E%s%s" % (printer._print(self.args[0]), printer._print(self.args[1]))
    
    def _pretty(self, printer, *args):
        return prettyForm(pretty_symbol("E_%s%s" % (printer._print(self.args[0]), printer._print(self.args[1]))))

class Corxet(sp.Expr):
    is_commutative = False
    def __new__(cls, A, B):
        i, j, k, l = A.i, A.j, B.i, B.j
        l = rl.relE[i][j][k][l]
        args = []
        for it in l:
            elem = sp.Mul(sp.S(it[0]), Eel(it[1], it[2]))
            args.append(elem)
        print(sp.srepr(sp.Add(*args)))
        return sp.Add(*args)
    
if __name__ == "__main__":
    E11=Eel(1, 1)
    E42=Eel(4, 2)
    E21=Eel(2, 1)
    print(Corxet(E11, E21))
    ec=Corxet(E21, E42)
    print(ec)
