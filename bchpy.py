#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 01-12-2020
# alex
# bchpy.py

import numpy as np
import sympy as sp
import relations as rl
import recursive as rc
import time as tm
import datetime
from sympy import re, im
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


class Metode():
    depth = 0
    w = []
    setw = False
    cofs = []
    
    def __init__(self, depth = 6):
        self.depth = len(rl.tamE) if depth >= len(rl.tamE) else depth
        self.w = self.__init_mat()
        
    def setBAB(self, *args, debug = False):
        if self.setw == True:
            self.w = self.__init_mat()
        self.cofs = list(args)
        cofsR = list(reversed(args))
        self.w[1][2] = cofsR[0]
        for i in range(1, len(cofsR)):
            if debug == True:
                txt_p = "Iteració " + str(i) + " amb Eel(1, " + str(2 - (i%2)) + ")"
                printd(txt_p)
            self.__recur_AB(cofsR[i], (i - 1)%2)
        self.setw = True

    def setABA(self, *args, debug = False):
        if self.setw == True:
            self.w = self.__init_mat()
        self.cofs = list(args)
        cofsR = list(reversed(args))
        self.w[1][1] = cofsR[0]
        for i in range(1, len(cofsR)):
            if debug == True:
                printd("Iteració " + str(i) + " amb Eel(1, " + str((i%2) + 1) + ")")
            self.__recur_AB(cofsR[i], i%2)
        self.setw = True

    def collectI(self):
        if self.setw == False:
            printe("Les equacions del mètode no estan calculades.")
            exit(-1)
        for i in range(2, len(self.w), 2):
            for j in range(1, len(self.w[i])):
                self.w[i][j] = self.w[i][j].collect(sp.I)

    def toExpr(self):
        esq = sp.S(0)
        for i in range(1, len(self.w)):
            for j in range(1, len(self.w[i])):
                esq += self.w[i][j]*Eel(i, j)
        return esq

    def importFromExpr(self, esq, debug = False):
        _x_diff_var = sp.Symbol("_x_diff_var")
        for i in range(self.depth):
            for j in range(rl.tamE[i]):
                if debug == True:
                    printd("Important Eel(" + str(i + 1) + "," + str(j + 1) + ")")
                self.w[i + 1][j + 1] = sp.expand(sp.diff(esq.subs(Eel(i + 1, j + 1), _x_diff_var), _x_diff_var)).subs(_x_diff_var, Eel(i + 1, j + 1))
        self.setw = True  

    def save(self, nom_fitxer):
        f = open(nom_fitxer, "w")
        f.write(str(self.depth) + "\n")
        f.write(str(self.cofs) + "\n")
        for i in range(self.depth + 1):
            f.write(str(self.w[i]))
            f.write("\n")
        f.close()

    def read(self, nom_fitxer):
        f = open(nom_fitxer, "r")
        linies = f.readlines()
        self.depth = int(linies[0].replace("\n", ""))
        cofs_txt = list(linies[1].replace("\n", "").replace("[", "").replace("]", "").split(","))
        self.cofs = []
        for i in range(len(cofs_txt)):
            self.cofs.append(sp.sympify(cofs_txt[i]))
        self.w = self.__init_mat()
        linies_w = linies[2:]
        for i in range(len(linies_w)):
            lin_txt = list(linies_w[i].replace("\n", "").replace("[", "").replace("]", "").split(","))
            for j in range(len(lin_txt)):
                self.w[i][j] = sp.sympify(lin_txt[j])
        f.close()
        
    def exportpy(self, nom_fitxer):
        if self.setw == False:
            printe("Les equacions del mètode no estan calculades.")
            exit(-1)        
        ara = datetime.datetime.now()
        f = open(nom_fitxer, "w")
        f.write("#!/usr/bin/env python\n# -*- coding: utf-8 -*-\n")
        f.write("# " + str(ara.day) + "-" + str(ara.month) + "-" + str(ara.year) + "\n")
        f.write("# " + nom_fitxer + "\n")
        f.write("# Fitxer generat automàticament per bchpy. No tocar!\n")
        cof_txt = str(self.cofs).replace("[", "").replace("]", "")
        f.write("# Esquema: " + cof_txt + "\n\n")
        for i in range(1, len(self.w)):
            for j in range(1, len(self.w[i])):
                cad = str(self.w[i][j])
                cad = strsub(r'a([0-9]*)', r'a[\1]', cad)
                cad = strsub(r'b([0-9]*)', r'br[\1]', cad)
                cad = strsub(r're\(br\[([0-9]*)\]\)', r'br[\1]', cad)
                cad = strsub(r'im\(br\[([0-9]*)\]\)', r'bi[\1]', cad)
                cad = cad.replace("I*", "")
                f.write("def w" + str(i) + str(j) + "(a, br, bi):\n")
                f.write("    return " + cad + "\n\n")
        f.close()
        
    def cprint(self):
        for i in range(1, len(self.w)):
            ctxt = "green"
            if i % 2 == 0:
                ctxt = "cyan"
            for j in range(1, len(self.w[i])):
                txt_p = "w(" + str(i) + "," + str(j) + ")"
                print(colored(txt_p, ctxt, attrs=['bold']), "=", self.w[i][j])

    def latex_str(self):
        cad = "\\begin{align}\n"
        for i in range(1, len(self.w)):
            for j in range(1, len(self.w[i])):
                cad = cad + "\omega_{" + str(i) + "," + str(j) + "} &= "+ sp.latex(self.w[i][j]) + "\\\ \n"
        cad = cad + "\end{align}"
        return cad

    def _latex(self, printer, *args):
        return self.latex_str()
        
    def __init_mat(self):
        wret = []
        for i in range(self.depth):
            wret.append([])
            for j in range(rl.tamE[i]):
                wret[i].append(sp.S(0))
            wret[i].insert(0, sp.S(0))
        wret.insert(0, [sp.S(0), sp.S(0)])
        return wret
    
    def __str__(self):
        cad = ""
        for i in range(1, len(self.w)):
            for j in range(1, len(self.w[i])):
                cad = cad + "w(" + str(i) + "," + str(j) + ") = " + str(self.w[i][j]) + "\n"
        return cad
    
    def __recur_AB(self, x, AB):
        bet = self.__init_mat()
        if AB == 1:
            bet = rc.recB(self.w, bet, x, self.depth)
        else:
            bet = rc.recA(self.w, bet, x, self.depth)               
        for i in range(len(bet)):
            for j in range(len(bet[i])):
                bet[i][j] = bet[i][j].expand()
        self.w = bet.copy()

def printd(text):
    print(colored("DEBUG:", "yellow", attrs=['bold']), text)
    
def printe(text):
    print(colored("ERROR:", "red", attrs=['bold']), text)

def bch9(A, B, depth = 6, debug = False):
    if debug == True:
        printd("BCH d'ordre 1 a 4")
    e21 = Corxet(A, B)
    e31 = Corxet(A, e21)
    e32 = Corxet(B, e21)
    e41 = Corxet(A, e31)
    e42 = Corxet(A, e32)
    e43 = Corxet(B, -e32)
    D = A + B
    D += sp.Rational(1, 2)*e21
    D += sp.Rational(1, 12)*(e31 - e32)
    D += sp.Rational(-1, 24)*e42
    if (depth >= 5):
        if debug == True:
            printd("BCH d'ordre 5")
        e51 = Corxet(A, e41)
        e52 = Corxet(B, e41)
        e53 = Corxet(A, -e42)
        e54 = Corxet(B, e42)
        e55 = Corxet(A, e43)
        e56 = Corxet(B, e43)
        D += sp.Rational(1, 720)*(-e51 - e56 + sp.S(6)*e53 + sp.S(6)*e54 + sp.S(2)*e55 + sp.S(2)*e52)
    if (depth >= 6):
        if debug == True:
            printd("BCH d'ordre 6")
        f6 = [sp.S(0)]*4
        f6[0] = Corxet(A, Corxet(A, Corxet(B, Corxet(B, Corxet(A, B))))) 
        f6[1] = Corxet(A, Corxet(B, Corxet(B, Corxet(A, Corxet(A, B)))))
        f6[2] = Corxet(A, Corxet(B, Corxet(B, Corxet(B, Corxet(A, B)))))
        f6[3] = Corxet(B, Corxet(A, Corxet(A, Corxet(A, Corxet(A, B)))))      
        den = sp.S([-720, 240, 1440, 1440])
        c6 = sp.S(0)
        for i in range(0, len(den)):
            c6 += sp.Rational(1, den[i])*f6[i]
        D += c6
    if (depth >= 7):
        if debug == True:
            printd("BCH d'ordre 7")
        f7 = [sp.S(0)]*18
        f7[0]  = Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(A, B))))))
        f7[1]  = Corxet(A, Corxet(A, Corxet(B, Corxet(A, Corxet(A, Corxet(A, B))))))
        f7[2]  = Corxet(A, Corxet(A, Corxet(B, Corxet(B, Corxet(A, Corxet(A, B))))))
        f7[3]  = Corxet(A, Corxet(B, Corxet(A, Corxet(A, Corxet(A, Corxet(A, B)))))) 
        f7[4]  = Corxet(A, Corxet(B, Corxet(A, Corxet(B, Corxet(A, Corxet(A, B))))))
        f7[5]  = Corxet(A, Corxet(B, Corxet(A, Corxet(B, Corxet(B, Corxet(A, B)))))) 
        f7[6]  = Corxet(A, Corxet(B, Corxet(B, Corxet(A, Corxet(A, Corxet(A, B))))))
        f7[7]  = Corxet(A, Corxet(B, Corxet(B, Corxet(B, Corxet(A, Corxet(A, B))))))
        f7[8]  = Corxet(A, Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(A, B))))))
        f7[9]  = Corxet(B, Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(A, B)))))) 
        f7[10] = Corxet(B, Corxet(A, Corxet(B, Corxet(A, Corxet(A, Corxet(A, B))))))
        f7[11] = Corxet(B, Corxet(A, Corxet(B, Corxet(B, Corxet(A, Corxet(A, B)))))) 
        f7[12] = Corxet(B, Corxet(B, Corxet(A, Corxet(A, Corxet(A, Corxet(A, B)))))) 
        f7[13] = Corxet(B, Corxet(B, Corxet(A, Corxet(B, Corxet(A, Corxet(A, B)))))) 
        f7[14] = Corxet(B, Corxet(B, Corxet(A, Corxet(B, Corxet(B, Corxet(A, B)))))) 
        f7[15] = Corxet(B, Corxet(B, Corxet(B, Corxet(A, Corxet(A, Corxet(A, B))))))
        f7[16] = Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(A, Corxet(A, B))))))
        f7[17] = Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(A, B))))))
        den = sp.S([30240, 5040, -10080, 10080, 1008, 5040, -7560, 3360, 10080, -10080, -1260, -1680, 3360, -3360, -2520, 7560, 10080, -30240])
        c7 = sp.S(0)
        for i in range(0, len(den)):
            c7 += sp.Rational(1, den[i])*f7[i]
        D += c7
    if (depth >= 8):
        if debug == True:
            printd("BCH d'ordre 8")
        f8 = [sp.S(0)]*13
        f8[0]  = Corxet(A, Corxet(A, Corxet(A, Corxet(B, Corxet(B, Corxet(B, Corxet(A, B)))))))
        f8[1]  = Corxet(A, Corxet(A, Corxet(B, Corxet(A, Corxet(B, Corxet(B, Corxet(A, B)))))))
        f8[2]  = Corxet(A, Corxet(A, Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(A, B)))))))
        f8[3]  = Corxet(A, Corxet(B, Corxet(A, Corxet(A, Corxet(B, Corxet(B, Corxet(A, B)))))))
        f8[4]  = Corxet(A, Corxet(B, Corxet(A, Corxet(B, Corxet(B, Corxet(A, Corxet(A, B))))))) 
        f8[5]  = Corxet(A, Corxet(B, Corxet(A, Corxet(B, Corxet(B, Corxet(B, Corxet(A, B)))))))
        f8[6]  = Corxet(A, Corxet(B, Corxet(B, Corxet(A, Corxet(B, Corxet(A, Corxet(A, B)))))))
        f8[7]  = Corxet(A, Corxet(B, Corxet(B, Corxet(A, Corxet(B, Corxet(B, Corxet(A, B)))))))
        f8[8]  = Corxet(A, Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(A, B)))))))
        f8[9]  = Corxet(B, Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(A, B)))))))
        f8[10] = Corxet(B, Corxet(A, Corxet(A, Corxet(A, Corxet(B, Corxet(A, Corxet(A, B)))))))
        f8[11] = Corxet(B, Corxet(A, Corxet(A, Corxet(B, Corxet(A, Corxet(A, Corxet(A, B)))))))
        f8[12] = Corxet(B, Corxet(B, Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(A, B)))))))
        den = sp.S([sp.Rational(-24192, 5), 2520, 20160, 15120, -2016, -20160, 20160, -10080, -60480, -60480, 20160, -5040, 20160])
        c8 = sp.S(0)
        for i in range(0, len(den)):
            c8 += sp.Rational(1, den[i])*f8[i]
        D += c8
    if (depth >= 9):
        if debug == True:
            printd("BCH d'ordre 9")
        f9 = [sp.S(0)]*38
        f9[0]  = Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(A, B))))))))
        f9[1]  = Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(B, Corxet(A, Corxet(A, B))))))))
        f9[2]  = Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(B, Corxet(B, Corxet(A, B))))))))
        f9[3]  = Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(B, Corxet(B, Corxet(A, Corxet(A, B))))))))
        f9[4]  = Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(B, Corxet(B, Corxet(B, Corxet(A, B)))))))) 
        f9[5]  = Corxet(A, Corxet(A, Corxet(A, Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(A, B))))))))
        f9[6]  = Corxet(A, Corxet(A, Corxet(B, Corxet(A, Corxet(A, Corxet(B, Corxet(A, Corxet(A, B))))))))
        f9[7]  = Corxet(A, Corxet(A, Corxet(B, Corxet(A, Corxet(A, Corxet(B, Corxet(B, Corxet(A, B))))))))
        f9[8]  = Corxet(A, Corxet(A, Corxet(B, Corxet(A, Corxet(B, Corxet(A, Corxet(A, Corxet(A, B))))))))
        f9[9]  = Corxet(A, Corxet(A, Corxet(B, Corxet(A, Corxet(B, Corxet(B, Corxet(B, Corxet(A, B))))))))
        f9[10] = Corxet(A, Corxet(A, Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(A, B))))))))
        f9[11] = Corxet(A, Corxet(B, Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(A, B))))))))
        f9[12] = Corxet(A, Corxet(B, Corxet(A, Corxet(A, Corxet(A, Corxet(B, Corxet(A, Corxet(A, B))))))))
        f9[13] = Corxet(A, Corxet(B, Corxet(A, Corxet(A, Corxet(B, Corxet(B, Corxet(A, Corxet(A, B))))))))
        f9[14] = Corxet(A, Corxet(B, Corxet(A, Corxet(A, Corxet(B, Corxet(B, Corxet(B, Corxet(A, B))))))))
        f9[15] = Corxet(A, Corxet(B, Corxet(A, Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(A, B))))))))
        f9[16] = Corxet(A, Corxet(B, Corxet(B, Corxet(A, Corxet(B, Corxet(B, Corxet(A, Corxet(A, B))))))))
        f9[17] = Corxet(A, Corxet(B, Corxet(B, Corxet(A, Corxet(B, Corxet(B, Corxet(B, Corxet(A, B))))))))
        f9[18] = Corxet(A, Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(A, Corxet(A, Corxet(A, B))))))))
        f9[19] = Corxet(A, Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(A, B))))))))
        f9[20] = Corxet(B, Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(A, B))))))))
        f9[21] = Corxet(B, Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(B, Corxet(A, Corxet(A, B))))))))
        f9[22] = Corxet(B, Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(B, Corxet(B, Corxet(A, B))))))))
        f9[23] = Corxet(B, Corxet(A, Corxet(A, Corxet(B, Corxet(A, Corxet(B, Corxet(A, Corxet(A, B))))))))
        f9[24] = Corxet(B, Corxet(A, Corxet(B, Corxet(B, Corxet(A, Corxet(A, Corxet(A, Corxet(A, B))))))))
        f9[25] = Corxet(B, Corxet(A, Corxet(B, Corxet(B, Corxet(A, Corxet(B, Corxet(A, Corxet(A, B))))))))
        f9[26] = Corxet(B, Corxet(A, Corxet(B, Corxet(B, Corxet(A, Corxet(B, Corxet(B, Corxet(A, B))))))))
        f9[27] = Corxet(B, Corxet(A, Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(A, B))))))))
        f9[28] = Corxet(B, Corxet(B, Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(A, B))))))))
        f9[29] = Corxet(B, Corxet(B, Corxet(A, Corxet(B, Corxet(A, Corxet(A, Corxet(A, Corxet(A, B))))))))
        f9[30] = Corxet(B, Corxet(B, Corxet(A, Corxet(B, Corxet(B, Corxet(A, Corxet(A, Corxet(A, B))))))))
        f9[31] = Corxet(B, Corxet(B, Corxet(A, Corxet(B, Corxet(B, Corxet(B, Corxet(A, Corxet(A, B))))))))
        f9[32] = Corxet(B, Corxet(B, Corxet(B, Corxet(A, Corxet(A, Corxet(A, Corxet(A, Corxet(A, B))))))))
        f9[33] = Corxet(B, Corxet(B, Corxet(B, Corxet(A, Corxet(A, Corxet(B, Corxet(B, Corxet(A, B))))))))
        f9[34] = Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(A, Corxet(A, Corxet(A, Corxet(A, B))))))))
        f9[35] = Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(A, Corxet(A, Corxet(A, B))))))))
        f9[36] = Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(A, Corxet(A, B))))))))
        f9[37] = Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(B, Corxet(A, B))))))))
        den = sp.S([-1209600, -302400, 113400, -40320, -120960, sp.Rational(604800, 13), 24192, 30240, -60480, -20160, sp.Rational(-604800, 11), -151200, -20160, -10080, 40320, 18900, -20160, -17280, -120960, -302400, 302400, 50400, 120960, 20160, -40320, 10080, 30240, 151200, 302400, 20160, -30240, 60480, sp.Rational(-604800, 13), -90720, 120960, 453600, 302400, 1209600])
        c9 = sp.S(0)
        for i in range(0, len(den)):
            c9 += sp.Rational(1, den[i])*f9[i]
        D += c9
    return D
