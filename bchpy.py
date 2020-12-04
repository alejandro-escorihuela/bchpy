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
    is_number = False
    is_zero = False
    def __init__(self, i = 1, j = 1):
        self.i = i
        self.j = j
        if i == 0:
            self.is_zero = True
            
    def _sympystr(self, printer, *args):
        if self.is_zero:
            return "0"
        return "E%s%s" % (printer._print(self.args[0]), printer._print(self.args[1]))
    
    def _pretty(self, printer, *args):
        if self.is_zero:
            return prettyForm(pretty_symbol("0"))
        return prettyForm(pretty_symbol("E_%s%s" % (printer._print(self.args[0]), printer._print(self.args[1]))))

class Corxet(sp.Expr):
    def __new__(cls, A, B):
        return cls.eval(cls, A, B)

    def eval(cls, A, B):
        if A.compare(B) == 1:
            return sp.S.NegativeOne*cls.eval(cls, B, A)
        i, j, k, l = A.i, A.j, B.i, B.j
        l = rl.corxetErel(i, j, k, l)
        args = []
        for it in l:
            elem = sp.Mul(sp.Rational(it[0], it[1]), Eel(it[2], it[3]))
            args.append(elem)
        return sp.Add(*args)

if __name__ == "__main__":
    print(Eel(0,0))
    E11=Eel(1, 1)
    E62=Eel(6, 2)
    E21=Eel(2, 1)
    print(Corxet(E11, E21))
    ec=Corxet(Eel(7,1), Eel(3,4))
    print(ec)
