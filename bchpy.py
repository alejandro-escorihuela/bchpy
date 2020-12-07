#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 01-12-2020
# alex
# bchpy.py

import numpy as np
import sympy as sp
import relations as rl
import time as tm
from sympy.printing.pretty.stringpict import prettyForm
from sympy.printing.pretty.pretty_symbology import pretty_symbol
from sympy.core.containers import Tuple

class Eel(sp.Expr):
    is_commutative = False
    is_number = False
    is_zero = False

    def __new__(cls, *args, **kwargs):
        args = cls._eval_args(args, **kwargs)
        return sp.Expr.__new__(cls, *args)
    
    @classmethod
    def _eval_args(cls, args):
        return tuple(Tuple(*args))
 
    def __init__(self, i = 1, j = 1):
        self.i = i
        self.j = j
        if i == 0 or j == 0:
            self.is_zero = True
    
    def _sympystr(self, printer, *args):
        if self.is_zero:
            return "0"
        return "E%s%s" % (printer._print(self.args[0]), printer._print(self.args[1]))
    
    def _sympyrepr(self, printer, *args):
        nom = self.__class__.__name__
        return "%s(%d,%d)" % (nom, self.i, self.j)
    
    def _pretty(self, printer, *args):
        if self.is_zero:
            return prettyForm(pretty_symbol("0"))
        return prettyForm(pretty_symbol("E_%s%s" % (printer._print(self.args[0]), printer._print(self.args[1]))))
    def _latex(self, printer, *args):
        return "E_{%d,%d}" % (self.i, self.j)
    
class Corxet(sp.Expr):
    def __new__(cls, A, B):
        return cls.eval(cls, A, B)

    def eval(cls, A, B):
        if A == B:
            return sp.S.Zero
        if A.is_commutative or B.is_commutative:
            return sp.S.Zero
        if isinstance(A, sp.Add):
            sargs = []
            for term in A.args:
                sargs.append(Corxet(term, B))
            return sp.Add(*sargs)
        if isinstance(B, sp.Add):
            sargs = []
            for term in B.args:
                sargs.append(Corxet(A, term))
            return sp.Add(*sargs)
        cA, ncA = A.args_cnc()
        cB, ncB = B.args_cnc()
        c_part = cA + cB
        if c_part:
            return sp.Mul(sp.Mul(*c_part), cls(sp.Mul._from_args(ncA), sp.Mul._from_args(ncB)))
        if A.compare(B) == 1:
            return sp.S.NegativeOne*cls.eval(cls, B, A)
        i, j, k, l = A.i, A.j, B.i, B.j
        l = rl.corxetErel(i, j, k, l)
        args = []
        for it in l:
            elem = sp.Mul(sp.Rational(it[0], it[1]), Eel(it[2], it[3]))
            args.append(elem)
        return sp.Add(*args)

def mat_sch(esq):
    w = []
    for i in range(len(rl.tamE)):
        w.append([])
        for j in range(rl.tamE[i]):
            # we = sp.simplify(sp.diff(esq.subs(Eel(i + 1, j + 1), x), x)).subs(x, Eel(i + 1, j + 1))
            we = sp.expand(sp.diff(esq.subs(Eel(i + 1, j + 1), x), x)).subs(x, Eel(i + 1, j + 1))
            w[i].append(we)
        w[i].insert(0, 0)
    w.insert(0, [0, 0])
    return w

def col_sch(esq):
    w = mat_sch(esq)
    esq_col = sp.S(0)
    for i in range(1, len(w)):
        for j in range(1, len(w[i])):
            esq_col += w[i][j]*Eel(i, j)
    return esq_col

def scheme(*args):
    cofs = list(reversed(args))
    esq = bch6(cofs[1]*Eel(1, 1), cofs[0]*Eel(1, 2))
    esq = col_sch(esq)
    for i in range(2, len(cofs)):
        print(i, 2 - (i%2))
        esq = bch6(cofs[i]*Eel(1, 2 - (i%2)), esq)
        esq = col_sch(esq)
    return esq

def bch6(A, B):
    e21 = Corxet(A, B)
    e31 = Corxet(A, e21)
    e32 = Corxet(B, e21)
    e41 = Corxet(A, e31)
    e42 = Corxet(A, e32)
    e43 = Corxet(B, -e32)
    e51 = Corxet(A, e41)
    e52 = Corxet(B, e41)
    e53 = Corxet(A, -e42)
    e54 = Corxet(B, e42)
    e55 = Corxet(A, e43)
    e56 = Corxet(B, e43)
    e62 = Corxet(B, e51)
    e64 = Corxet(A, e54)
    e66 = Corxet(A, e55)
    e68 = Corxet(A, e56)
    D = A + B
    D += sp.Rational(1, 2)*e21
    D += sp.Rational(1, 12)*(e31 - e32)
    D += sp.Rational(-1, 24)*e42
    D += sp.Rational(1, 720)*(-e51 - e56 + sp.S(6)*e53 + sp.S(6)*e54 + sp.S(2)*e55 + sp.S(2)*e52)
    D += sp.Rational(1, 1440)*(e62 - e68 + sp.S(2)*e66 + sp.S(6)*e64)
    return D

if __name__ == "__main__":
    print("BCH")
    x, h = sp.symbols("x h")
    a = sp.symbols("a0:3")
    b = sp.symbols("b0:4")
    t0 = tm.time()
    met = scheme(b[0], a[0], b[1], a[1], b[2], sp.Rational(1, 2) - a[1] - a[0], sp.S(1) - (sp.S(2)*(b[0] + b[1] + b[2])), sp.Rational(1, 2) - a[1] - a[0], b[2], a[1], b[1], a[0], b[0])
    w = mat_sch(met)
    for i in range(1, len(w)):
        for j in range(1, len(w[i])):
            print("w(" + str(i) + "," + str(j) + ") =", w[i][j])
    t1 = tm.time()
    print(t1-t0, "s")
