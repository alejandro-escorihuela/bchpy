#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 01-12-2020
# alex
# bchpy.py

import numpy as np
import sympy as sp
import relations as rl
import time as tm
import datetime
from sympy.printing.pretty.stringpict import prettyForm
from sympy.printing.pretty.pretty_symbology import pretty_symbol
from sympy.core.containers import Tuple
from sympy.parsing.sympy_parser import parse_expr
from termcolor import colored
from re import sub as strsub

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


def init_mat():
    w = []
    for i in range(len(rl.tamE)):
        w.append([])
        for j in range(rl.tamE[i]):
            w[i].append(sp.S(0))
        w[i].insert(0, sp.S(0))
    w.insert(0, [sp.S(0), sp.S(0)])
    return w

def esq2mat(esq):
    w = []
    _x_diff_var = sp.Symbol("_x_diff_var")
    for i in range(len(rl.tamE)):
        w.append([])
        for j in range(rl.tamE[i]):
            we = sp.expand(sp.diff(esq.subs(Eel(i + 1, j + 1), _x_diff_var), _x_diff_var)).subs(_x_diff_var, Eel(i + 1, j + 1))
            w[i].append(we)
        w[i].insert(0, 0)
    w.insert(0, [0, 0])
    return w

def mat2esq(w):
    esq = sp.S(0)
    for i in range(1, len(w)):
        for j in range(1, len(w[i])):
            esq += w[i][j]*Eel(i, j)
    return esq

def col_esq(esq):
    w = esq2mat(esq)
    esq_col = sp.S(0)
    for i in range(1, len(w)):
        for j in range(1, len(w[i])):
            esq_col += w[i][j]*Eel(i, j)
    return esq_col

def esqBAB(*args, debug = False):
    cofs = list(reversed(args))
    if debug == True:
        printd("Iteració 1 BCH(Eel(1, 1), Eel(1, 2))")
    esq = bch6(cofs[1]*Eel(1, 1), cofs[0]*Eel(1, 2))
    esq = col_esq(esq)
    for i in range(2, len(cofs)):
        if debug == True:
            txt_p = "Iteració " + str(i) + " amb Eel(1, " + str(2 - (i%2)) + ")"
            printd(txt_p)
        esq = bch6(cofs[i]*Eel(1, 2 - (i%2)), esq)
        esq = col_esq(esq)
    return esq

def mat_esqBAB(*args, debug = False):
    cofs = list(reversed(args))
    w = init_mat()
    w[1][2] = cofs[0]
    for i in range(1, len(cofs)):
        if debug == True:
            txt_p = "Iteració " + str(i) + " amb Eel(1, " + str(2 - (i%2)) + ")"
            printd(txt_p)
        w = recur_AB(cofs[i], w, (i - 1)%2)
    return w

def mat_esqABA(*args, debug = False):
    cofs = list(reversed(args))
    w = init_mat()
    w[1][1] = cofs[0]
    for i in range(1, len(cofs)):
        if debug == True:
            printd("Iteració " + str(i) + " amb Eel(1, " + str((i%2) + 1) + ")")
        w = recur_AB(cofs[i], w, i%2)
    return w

def printd(text):
    print(colored("DEBUG:", "yellow", attrs=['bold']), text)

def print_mat(w):
    for i in range(1, len(w)):
        ctxt = "green"
        if i % 2 == 0:
            ctxt = "cyan"
        for j in range(1, len(w[i])):
            txt_p = "w(" + str(i) + "," + str(j) + ")"
            print(colored(txt_p, ctxt, attrs=['bold']), "=", w[i][j])
            
def mat_collectI(w):
    for i in range(2, len(w), 2):
        for j in range(1, len(w[i])):
            w[i][j] = w[i][j].collect(sp.I)
    return w

def mat_exportC(w, nom_fitxer):
    ara = datetime.datetime.now()
    f = open(nom_fitxer, "w")
    f.write("#!/usr/bin/env python\n# -*- coding: utf-8 -*-\n")
    f.write("# " + str(ara.day) + "-" + str(ara.month) + "-" + str(ara.year) + "\n")
    f.write("# Fitxer generat automàticament per bchpy. No tocar!\n")
    f.write("# " + nom_fitxer + "\n\n")
    for i in range(1, len(w)):
        for j in range(1, len(w[i])):
            cad = str(w[i][j])
            cad = strsub(r'a([0-9]*)', r'a[\1]', cad)
            cad = strsub(r're\(b([0-9]*)\)', r'br[\1]', cad)
            cad = strsub(r'im\(b([0-9]*)\)', r'bi[\1]', cad)
            cad = cad.replace("I*", "")
            f.write("def w" + str(i) + str(j) + "(a, br, bi):\n")
            f.write("    return " + cad + "\n\n")
    f.close()

def recur_AB(x, alp, AB):
    bet = init_mat()
    if AB == 1:
        bet[1][1] = alp[1][1]
        bet[1][2] = x + alp[1][2]
        bet[2][1] = -x*alp[1][1]/2 + alp[2][1]
        bet[3][1] = x*alp[1][1]**2/12 + alp[3][1]
        bet[3][2] = -x**2*alp[1][1]/12 + x*alp[1][1]*alp[1][2]/12 + x*alp[2][1]/2 + alp[3][2]
        bet[4][1] = alp[4][1]
        bet[4][2] = x**2*alp[1][1]**2/24 - x*alp[1][1]*alp[2][1]/12 + x*alp[3][1]/2 + alp[4][2]
        bet[4][3] = -x**2*alp[1][1]*alp[1][2]/24 - x**2*alp[2][1]/12 + x*alp[1][2]*alp[2][1]/12 - x*alp[3][2]/2 + alp[4][3]
        bet[5][1] = -x*alp[1][1]**4/720 + alp[5][1]
        bet[5][2] = x**2*alp[1][1]**3/360 - x*alp[1][1]**3*alp[1][2]/720 + x*alp[1][1]*alp[3][1]/12 + x*alp[4][1]/2 + alp[5][2]
        bet[5][3] = x**2*alp[1][1]**3/120 + x*alp[1][1]**3*alp[1][2]/360 + x*alp[1][1]*alp[3][1]/6 + alp[5][3]
        bet[5][4] = x**3*alp[1][1]**2/120 - x**2*alp[1][1]**2*alp[1][2]/360 - x**2*alp[1][1]*alp[2][1]/24 + x**2*alp[3][1]/12 - x*alp[1][1]**2*alp[1][2]**2/360 + x*alp[1][1]*alp[3][2]/12 - x*alp[1][2]*alp[3][1]/12 + x*alp[2][1]**2/12 + x*alp[4][2]/2 + alp[5][4]
        bet[5][5] = x**3*alp[1][1]**2/360 + x**2*alp[1][1]**2*alp[1][2]/120 + x*alp[1][1]**2*alp[1][2]**2/720 + x*alp[1][1]*alp[3][2]/6 + x*alp[2][1]**2/12 + alp[5][5]
        bet[5][6] = -x**4*alp[1][1]/720 - x**3*alp[1][1]*alp[1][2]/180 + x**2*alp[1][1]*alp[1][2]**2/180 + x**2*alp[1][2]*alp[2][1]/24 - x**2*alp[3][2]/12 + x*alp[1][1]*alp[1][2]**3/720 + x*alp[1][2]*alp[3][2]/12 + x*alp[4][3]/2 + alp[5][6]
        bet[6][1] = alp[6][1]
        bet[6][2] = -x**2*alp[1][1]**4/1440 + x*alp[1][1]**3*alp[2][1]/720 + x*alp[1][1]*alp[4][1]/12 + x*alp[5][1]/2 + alp[6][2]
        bet[6][3] = -x*alp[1][1]*alp[4][1]/6 + alp[6][3]
        bet[6][4] = -x**3*alp[1][1]**3/240 - x**2*alp[1][1]**3*alp[1][2]/720 + x**2*alp[1][1]**2*alp[2][1]/120 - x**2*alp[1][1]*alp[3][1]/12 + x*alp[1][1]**2*alp[1][2]*alp[2][1]/360 - x*alp[1][1]*alp[4][2]/12 + x*alp[2][1]*alp[3][1]/12 - x*alp[5][3]/2 + alp[6][4]
        bet[6][5] = -x**2*alp[1][1]**3*alp[1][2]/864 + x**2*alp[1][1]*alp[3][1]/72 + x**2*alp[4][1]/12 + x*alp[1][1]**2*alp[1][2]*alp[2][1]/432 + x*alp[1][1]*alp[4][2]/36 - x*alp[1][2]*alp[4][1]/12 + x*alp[5][2]/2 - x*alp[5][3]/6 + alp[6][5]
        bet[6][6] = -x**3*alp[1][1]**3/720 - x**2*alp[1][1]**3*alp[1][2]/2160 + x**2*alp[1][1]**2*alp[2][1]/360 - x**2*alp[1][1]*alp[3][1]/36 + x*alp[1][1]**2*alp[1][2]*alp[2][1]/1080 + x*alp[1][1]*alp[4][2]/36 + x*alp[2][1]*alp[3][1]/12 - x*alp[5][3]/6 + alp[6][6]
        bet[6][7] = x**3*alp[1][1]**2*alp[1][2]/144 + x**3*alp[1][1]*alp[2][1]/72 + x**2*alp[1][1]**2*alp[1][2]**2/288 - x**2*alp[1][1]*alp[1][2]*alp[2][1]/72 + x**2*alp[1][2]*alp[3][1]/12 - x**2*alp[2][1]**2/24 - x**2*alp[4][2]/6 - x*alp[1][1]*alp[1][2]**2*alp[2][1]/144 + x*alp[1][1]*alp[4][3]/12 + x*alp[1][2]*alp[4][2]/6 - x*alp[2][1]*alp[3][2]/12 - x*alp[5][4] + x*alp[5][5]/2 + alp[6][7]
        bet[6][8] = x**4*alp[1][1]**2/1440 - x**3*alp[1][1]**2*alp[1][2]/720 - x**3*alp[1][1]*alp[2][1]/120 - x**2*alp[1][1]**2*alp[1][2]**2/720 + x**2*alp[1][1]*alp[1][2]*alp[2][1]/360 + x**2*alp[1][1]*alp[3][2]/24 - x**2*alp[1][2]*alp[3][1]/24 + x**2*alp[2][1]**2/24 + x**2*alp[4][2]/12 + x*alp[1][1]*alp[1][2]**2*alp[2][1]/360 - x*alp[1][1]*alp[4][3]/6 - x*alp[1][2]*alp[4][2]/12 + x*alp[2][1]*alp[3][2]/12 + x*alp[5][4]/2 + alp[6][8]
        bet[6][9] = x**4*alp[1][1]*alp[1][2]/1440 + x**4*alp[2][1]/720 + x**3*alp[1][1]*alp[1][2]**2/360 + x**3*alp[1][2]*alp[2][1]/180 + x**2*alp[1][1]*alp[1][2]**3/1440 - x**2*alp[1][2]**2*alp[2][1]/180 + x**2*alp[1][2]*alp[3][2]/24 + x**2*alp[4][3]/12 - x*alp[1][2]**3*alp[2][1]/720 - x*alp[1][2]*alp[4][3]/12 + x*alp[5][6]/2 + alp[6][9]
    else:
        bet[1][1] = x + alp[1][1]
        bet[1][2] = alp[1][2]
        bet[2][1] = x*alp[1][2]/2 + alp[2][1]
        bet[3][1] = x**2*alp[1][2]/12 - x*alp[1][1]*alp[1][2]/12 + x*alp[2][1]/2 + alp[3][1]
        bet[3][2] = -x*alp[1][2]**2/12 + alp[3][2]
        bet[4][1] = -x**2*alp[1][1]*alp[1][2]/24 + x**2*alp[2][1]/12 - x*alp[1][1]*alp[2][1]/12 + x*alp[3][1]/2 + alp[4][1]
        bet[4][2] = -x**2*alp[1][2]**2/24 - x*alp[1][2]*alp[2][1]/12 + x*alp[3][2]/2 + alp[4][2]
        bet[4][3] = alp[4][3]
        bet[5][1] = -x**4*alp[1][2]/720 - x**3*alp[1][1]*alp[1][2]/180 + x**2*alp[1][1]**2*alp[1][2]/180 - x**2*alp[1][1]*alp[2][1]/24 + x**2*alp[3][1]/12 + x*alp[1][1]**3*alp[1][2]/720 - x*alp[1][1]*alp[3][1]/12 + x*alp[4][1]/2 + alp[5][1]
        bet[5][2] = x**3*alp[1][2]**2/360 + x**2*alp[1][1]*alp[1][2]**2/120 + x*alp[1][1]**2*alp[1][2]**2/720 - x*alp[1][2]*alp[3][1]/6 + x*alp[2][1]**2/12 + alp[5][2]
        bet[5][3] = x**3*alp[1][2]**2/120 - x**2*alp[1][1]*alp[1][2]**2/360 + x**2*alp[1][2]*alp[2][1]/24 - x**2*alp[3][2]/12 - x*alp[1][1]**2*alp[1][2]**2/360 + x*alp[1][1]*alp[3][2]/12 - x*alp[1][2]*alp[3][1]/12 + x*alp[2][1]**2/12 - x*alp[4][2]/2 + alp[5][3]
        bet[5][4] = x**2*alp[1][2]**3/120 + x*alp[1][1]*alp[1][2]**3/360 - x*alp[1][2]*alp[3][2]/6 + alp[5][4]
        bet[5][5] = x**2*alp[1][2]**3/360 - x*alp[1][1]*alp[1][2]**3/720 - x*alp[1][2]*alp[3][2]/12 + x*alp[4][3]/2 + alp[5][5]
        bet[5][6] = -x*alp[1][2]**4/720 + alp[5][6]
        bet[6][1] = x**4*alp[1][1]*alp[1][2]/1440 - x**4*alp[2][1]/720 + x**3*alp[1][1]**2*alp[1][2]/360 - x**3*alp[1][1]*alp[2][1]/180 + x**2*alp[1][1]**3*alp[1][2]/1440 + x**2*alp[1][1]**2*alp[2][1]/180 - x**2*alp[1][1]*alp[3][1]/24 + x**2*alp[4][1]/12 + x*alp[1][1]**3*alp[2][1]/720 - x*alp[1][1]*alp[4][1]/12 + x*alp[5][1]/2 + alp[6][1]
        bet[6][2] = x**4*alp[1][2]**2/1440 - x**3*alp[1][1]*alp[1][2]**2/720 + x**3*alp[1][2]*alp[2][1]/120 - x**2*alp[1][1]**2*alp[1][2]**2/720 - x**2*alp[1][1]*alp[1][2]*alp[2][1]/360 + x**2*alp[1][1]*alp[3][2]/24 - x**2*alp[1][2]*alp[3][1]/24 + x**2*alp[2][1]**2/24 - x**2*alp[4][2]/12 - x*alp[1][1]**2*alp[1][2]*alp[2][1]/360 + x*alp[1][1]*alp[4][2]/12 - x*alp[1][2]*alp[4][1]/6 + x*alp[2][1]*alp[3][1]/12 + x*alp[5][3]/2 + alp[6][2]
        bet[6][3] = x**3*alp[1][1]*alp[1][2]**2/144 - x**3*alp[1][2]*alp[2][1]/72 + x**2*alp[1][1]**2*alp[1][2]**2/288 + x**2*alp[1][1]*alp[1][2]*alp[2][1]/72 - x**2*alp[1][1]*alp[3][2]/12 - x**2*alp[2][1]**2/24 + x**2*alp[4][2]/6 + x*alp[1][1]**2*alp[1][2]*alp[2][1]/144 - x*alp[1][1]*alp[4][2]/6 + x*alp[1][2]*alp[4][1]/12 - x*alp[2][1]*alp[3][1]/12 + x*alp[5][2]/2 - x*alp[5][3] + alp[6][3]
        bet[6][4] = x**3*alp[1][2]**3/240 + x**2*alp[1][1]*alp[1][2]**3/720 + x**2*alp[1][2]**2*alp[2][1]/120 - x**2*alp[1][2]*alp[3][2]/12 + x*alp[1][1]*alp[1][2]**2*alp[2][1]/360 - x*alp[1][2]*alp[4][2]/12 - x*alp[2][1]*alp[3][2]/12 + x*alp[5][4]/2 + alp[6][4]
        bet[6][5] = -x*alp[1][2]*alp[4][2]/18 + x*alp[2][1]*alp[3][2]/18 + alp[6][5]
        bet[6][6] = x**3*alp[1][2]**3/720 - x**2*alp[1][1]*alp[1][2]**3/1440 + x**2*alp[1][2]**2*alp[2][1]/360 - x**2*alp[1][2]*alp[3][2]/24 + x**2*alp[4][3]/12 - x*alp[1][1]*alp[1][2]**2*alp[2][1]/720 - x*alp[1][1]*alp[4][3]/12 - x*alp[1][2]*alp[4][2]/18 - x*alp[2][1]*alp[3][2]/36 + x*alp[5][5]/2 + alp[6][6]
        bet[6][7] = -x*alp[1][2]*alp[4][3]/6 + alp[6][7]
        bet[6][8] = -x**2*alp[1][2]**4/1440 - x*alp[1][2]**3*alp[2][1]/720 + x*alp[1][2]*alp[4][3]/12 + x*alp[5][6]/2 + alp[6][8]
        bet[6][9] = alp[6][9]
    for i in range(len(bet)):
        for j in range(len(bet[i])):
            bet[i][j] = bet[i][j].expand()
    return bet;

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
