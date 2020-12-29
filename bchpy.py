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

    def importFromExpr(self, esq):
        _x_diff_var = sp.Symbol("_x_diff_var")
        for i in range(self.depth):
            for j in range(rl.tamE[i]):
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
                
    def latex(self):
        print(self.latex_str())
        
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
            bet[1][1] = self.w[1][1]
            bet[1][2] = x + self.w[1][2]
            bet[2][1] = -x*self.w[1][1]/2 + self.w[2][1]
            bet[3][1] = x*self.w[1][1]**2/12 + self.w[3][1]
            bet[3][2] = -x**2*self.w[1][1]/12 + x*self.w[1][1]*self.w[1][2]/12 + x*self.w[2][1]/2 + self.w[3][2]
            bet[4][1] = self.w[4][1]
            bet[4][2] = x**2*self.w[1][1]**2/24 - x*self.w[1][1]*self.w[2][1]/12 + x*self.w[3][1]/2 + self.w[4][2]
            bet[4][3] = -x**2*self.w[1][1]*self.w[1][2]/24 - x**2*self.w[2][1]/12 + x*self.w[1][2]*self.w[2][1]/12 - x*self.w[3][2]/2 + self.w[4][3]
            if self.depth >= 5:
                bet[5][1] = -x*self.w[1][1]**4/720 + self.w[5][1]
                bet[5][2] = x**2*self.w[1][1]**3/360 - x*self.w[1][1]**3*self.w[1][2]/720 + x*self.w[1][1]*self.w[3][1]/12 + x*self.w[4][1]/2 + self.w[5][2]
                bet[5][3] = x**2*self.w[1][1]**3/120 + x*self.w[1][1]**3*self.w[1][2]/360 + x*self.w[1][1]*self.w[3][1]/6 + self.w[5][3]
                bet[5][4] = x**3*self.w[1][1]**2/120 - x**2*self.w[1][1]**2*self.w[1][2]/360 - x**2*self.w[1][1]*self.w[2][1]/24 + x**2*self.w[3][1]/12 - x*self.w[1][1]**2*self.w[1][2]**2/360 + x*self.w[1][1]*self.w[3][2]/12 - x*self.w[1][2]*self.w[3][1]/12 + x*self.w[2][1]**2/12 + x*self.w[4][2]/2 + self.w[5][4]
                bet[5][5] = x**3*self.w[1][1]**2/360 + x**2*self.w[1][1]**2*self.w[1][2]/120 + x*self.w[1][1]**2*self.w[1][2]**2/720 + x*self.w[1][1]*self.w[3][2]/6 + x*self.w[2][1]**2/12 + self.w[5][5]
                bet[5][6] = -x**4*self.w[1][1]/720 - x**3*self.w[1][1]*self.w[1][2]/180 + x**2*self.w[1][1]*self.w[1][2]**2/180 + x**2*self.w[1][2]*self.w[2][1]/24 - x**2*self.w[3][2]/12 + x*self.w[1][1]*self.w[1][2]**3/720 + x*self.w[1][2]*self.w[3][2]/12 + x*self.w[4][3]/2 + self.w[5][6]
            if self.depth >= 6:
                bet[6][1] = self.w[6][1]
                bet[6][2] = -x**2*self.w[1][1]**4/1440 + x*self.w[1][1]**3*self.w[2][1]/720 + x*self.w[1][1]*self.w[4][1]/12 + x*self.w[5][1]/2 + self.w[6][2]
                bet[6][3] = -x*self.w[1][1]*self.w[4][1]/6 + self.w[6][3]
                bet[6][4] = -x**3*self.w[1][1]**3/240 - x**2*self.w[1][1]**3*self.w[1][2]/720 + x**2*self.w[1][1]**2*self.w[2][1]/120 - x**2*self.w[1][1]*self.w[3][1]/12 + x*self.w[1][1]**2*self.w[1][2]*self.w[2][1]/360 - x*self.w[1][1]*self.w[4][2]/12 + x*self.w[2][1]*self.w[3][1]/12 - x*self.w[5][3]/2 + self.w[6][4]
                bet[6][5] = -x**2*self.w[1][1]**3*self.w[1][2]/864 + x**2*self.w[1][1]*self.w[3][1]/72 + x**2*self.w[4][1]/12 + x*self.w[1][1]**2*self.w[1][2]*self.w[2][1]/432 + x*self.w[1][1]*self.w[4][2]/36 - x*self.w[1][2]*self.w[4][1]/12 + x*self.w[5][2]/2 - x*self.w[5][3]/6 + self.w[6][5]
                bet[6][6] = -x**3*self.w[1][1]**3/720 - x**2*self.w[1][1]**3*self.w[1][2]/2160 + x**2*self.w[1][1]**2*self.w[2][1]/360 - x**2*self.w[1][1]*self.w[3][1]/36 + x*self.w[1][1]**2*self.w[1][2]*self.w[2][1]/1080 + x*self.w[1][1]*self.w[4][2]/36 + x*self.w[2][1]*self.w[3][1]/12 - x*self.w[5][3]/6 + self.w[6][6]
                bet[6][7] = x**3*self.w[1][1]**2*self.w[1][2]/144 + x**3*self.w[1][1]*self.w[2][1]/72 + x**2*self.w[1][1]**2*self.w[1][2]**2/288 - x**2*self.w[1][1]*self.w[1][2]*self.w[2][1]/72 + x**2*self.w[1][2]*self.w[3][1]/12 - x**2*self.w[2][1]**2/24 - x**2*self.w[4][2]/6 - x*self.w[1][1]*self.w[1][2]**2*self.w[2][1]/144 + x*self.w[1][1]*self.w[4][3]/12 + x*self.w[1][2]*self.w[4][2]/6 - x*self.w[2][1]*self.w[3][2]/12 - x*self.w[5][4] + x*self.w[5][5]/2 + self.w[6][7]
                bet[6][8] = x**4*self.w[1][1]**2/1440 - x**3*self.w[1][1]**2*self.w[1][2]/720 - x**3*self.w[1][1]*self.w[2][1]/120 - x**2*self.w[1][1]**2*self.w[1][2]**2/720 + x**2*self.w[1][1]*self.w[1][2]*self.w[2][1]/360 + x**2*self.w[1][1]*self.w[3][2]/24 - x**2*self.w[1][2]*self.w[3][1]/24 + x**2*self.w[2][1]**2/24 + x**2*self.w[4][2]/12 + x*self.w[1][1]*self.w[1][2]**2*self.w[2][1]/360 - x*self.w[1][1]*self.w[4][3]/6 - x*self.w[1][2]*self.w[4][2]/12 + x*self.w[2][1]*self.w[3][2]/12 + x*self.w[5][4]/2 + self.w[6][8]
                bet[6][9] = x**4*self.w[1][1]*self.w[1][2]/1440 + x**4*self.w[2][1]/720 + x**3*self.w[1][1]*self.w[1][2]**2/360 + x**3*self.w[1][2]*self.w[2][1]/180 + x**2*self.w[1][1]*self.w[1][2]**3/1440 - x**2*self.w[1][2]**2*self.w[2][1]/180 + x**2*self.w[1][2]*self.w[3][2]/24 + x**2*self.w[4][3]/12 - x*self.w[1][2]**3*self.w[2][1]/720 - x*self.w[1][2]*self.w[4][3]/12 + x*self.w[5][6]/2 + self.w[6][9]
            if self.depth >= 7:
                bet[7][1] = x*self.w[1][1]**6/30240 + self.w[7][1]
                bet[7][2] = -x**2*self.w[1][1]**5/10080 + x*self.w[1][1]**5*self.w[1][2]/30240 - x*self.w[1][1]**3*self.w[3][1]/720 + x*self.w[1][1]*self.w[5][1]/12 + x*self.w[6][1]/2 + self.w[7][2]
                bet[7][3] = x**2*self.w[1][1]**5/10080 - x*self.w[1][1]**5*self.w[1][2]/30240 - x*self.w[1][1]*self.w[5][1]/6 + self.w[7][3]
                bet[7][4] = x**3*self.w[1][1]**4/3360 + x**2*self.w[1][1]**4*self.w[1][2]/5040 + x**2*self.w[1][1]**3*self.w[2][1]/1440 + x**2*self.w[1][1]**2*self.w[3][1]/72 + x**2*self.w[1][1]*self.w[4][1]/24 + x**2*self.w[5][1]/12 - x*self.w[1][1]**4*self.w[1][2]**2/30240 + x*self.w[1][1]**3*self.w[3][2]/720 + x*self.w[1][1]*self.w[5][3]/12 - x*self.w[1][2]*self.w[5][1]/12 - x*self.w[2][1]*self.w[4][1]/12 + x*self.w[3][1]**2/12 + x*self.w[6][2]/2 + self.w[7][4]
                bet[7][5] = x**2*self.w[1][1]**5/5040 + x*self.w[1][1]**5*self.w[1][2]/6048 + x*self.w[1][1]**3*self.w[3][1]/360 + self.w[7][5]
                bet[7][6] = -x**3*self.w[1][1]**4/1260 - x**2*self.w[1][1]**4*self.w[1][2]/3360 - x**2*self.w[1][1]**2*self.w[3][1]/40 - x**2*self.w[1][1]*self.w[4][1]/12 + x*self.w[1][1]**4*self.w[1][2]**2/6048 - x*self.w[1][1]**3*self.w[3][2]/360 - x*self.w[1][1]**2*self.w[1][2]*self.w[3][1]/720 - x*self.w[1][1]**2*self.w[2][1]**2/720 + x*self.w[1][1]*self.w[5][2]/12 - x*self.w[1][1]*self.w[5][3]/6 + x*self.w[2][1]*self.w[4][1]/4 - x*self.w[3][1]**2/6 + x*self.w[6][3]/2 + self.w[7][6]
                bet[7][7] = x**3*self.w[1][1]**4/1120 + x**2*self.w[1][1]**4*self.w[1][2]/1680 + x**2*self.w[1][1]**2*self.w[3][1]/40 + x*self.w[1][1]**4*self.w[1][2]**2/7560 - x*self.w[1][1]**3*self.w[3][2]/360 + x*self.w[1][1]**2*self.w[1][2]*self.w[3][1]/120 - x*self.w[1][1]**2*self.w[2][1]**2/720 + x*self.w[1][1]*self.w[5][3]/6 + x*self.w[3][1]**2/12 + self.w[7][7]
                bet[7][8] = -x**4*self.w[1][1]**3/1120 - x**3*self.w[1][1]**3*self.w[1][2]/10080 + x**3*self.w[1][1]**2*self.w[2][1]/240 - x**3*self.w[1][1]*self.w[3][1]/60 + x**2*self.w[1][1]**3*self.w[1][2]**2/5040 + x**2*self.w[1][1]**2*self.w[1][2]*self.w[2][1]/720 - x**2*self.w[1][1]**2*self.w[3][2]/360 + x**2*self.w[1][1]*self.w[1][2]*self.w[3][1]/180 - x**2*self.w[1][1]*self.w[2][1]**2/180 - x**2*self.w[1][1]*self.w[4][2]/24 + x**2*self.w[2][1]*self.w[3][1]/24 - x**2*self.w[5][3]/12 + x*self.w[1][1]**3*self.w[1][2]**3/7560 - x*self.w[1][1]**2*self.w[1][2]*self.w[3][2]/180 + x*self.w[1][1]*self.w[1][2]**2*self.w[3][1]/180 - x*self.w[1][1]*self.w[1][2]*self.w[2][1]**2/240 + x*self.w[1][1]*self.w[5][4]/12 + x*self.w[1][2]*self.w[5][3]/12 + x*self.w[2][1]*self.w[4][2]/6 - x*self.w[3][1]*self.w[3][2]/12 + x*self.w[6][4]/2 + self.w[7][8]
                bet[7][9] = x**3*self.w[1][1]**4/5040 + 11*x**2*self.w[1][1]**4*self.w[1][2]/30240 + x**2*self.w[1][1]**2*self.w[3][1]/360 + x*self.w[1][1]**4*self.w[1][2]**2/18144 + x*self.w[1][1]**3*self.w[3][2]/540 + x*self.w[1][1]**2*self.w[1][2]*self.w[3][1]/1080 + x*self.w[1][1]**2*self.w[2][1]**2/1080 - x*self.w[1][1]*self.w[5][2]/6 + x*self.w[1][1]*self.w[5][3]/18 - x*self.w[2][1]*self.w[4][1]/6 + x*self.w[3][1]**2/18 + self.w[7][9]
                bet[7][10] = x**4*self.w[1][1]**3/30240 - x**3*self.w[1][1]**3*self.w[1][2]/10080 + x**3*self.w[1][1]*self.w[3][1]/360 + x**2*self.w[1][1]**3*self.w[1][2]**2/5040 + x**2*self.w[1][1]**2*self.w[1][2]*self.w[2][1]/864 + x**2*self.w[1][1]**2*self.w[3][2]/540 - x**2*self.w[1][1]*self.w[1][2]*self.w[3][1]/1080 + x**2*self.w[1][1]*self.w[2][1]**2/1080 + x**2*self.w[1][1]*self.w[4][2]/72 - x**2*self.w[1][2]*self.w[4][1]/24 + x**2*self.w[5][2]/12 - x**2*self.w[5][3]/36 + x*self.w[1][1]**3*self.w[1][2]**3/18144 + x*self.w[1][1]**2*self.w[1][2]*self.w[3][2]/720 - x*self.w[1][1]*self.w[1][2]**2*self.w[3][1]/1080 - x*self.w[1][1]*self.w[1][2]*self.w[2][1]**2/2160 - x*self.w[1][2]*self.w[5][2]/12 + x*self.w[1][2]*self.w[5][3]/36 - x*self.w[2][1]*self.w[4][2]/36 + x*self.w[3][1]*self.w[3][2]/36 + x*self.w[6][5]/2 + self.w[7][10]
                bet[7][11] = x**3*self.w[1][1]**4/3024 + x**2*self.w[1][1]**4*self.w[1][2]/15120 + x**2*self.w[1][1]**2*self.w[3][1]/90 - x*self.w[1][1]**4*self.w[1][2]**2/90720 - x*self.w[1][1]**3*self.w[3][2]/270 + x*self.w[1][1]**2*self.w[1][2]*self.w[3][1]/270 - x*self.w[1][1]**2*self.w[2][1]**2/540 + x*self.w[1][1]*self.w[5][3]/18 + x*self.w[3][1]**2/18 + self.w[7][11]
                bet[7][12] = -x**4*self.w[1][1]**3/10080 + x**3*self.w[1][1]**3*self.w[1][2]/1890 + x**3*self.w[1][1]**2*self.w[2][1]/720 - x**3*self.w[1][1]*self.w[3][1]/360 + x**2*self.w[1][1]**3*self.w[1][2]**2/10080 + x**2*self.w[1][1]**2*self.w[1][2]*self.w[2][1]/2160 + 11*x**2*self.w[1][1]**2*self.w[3][2]/1080 + 11*x**2*self.w[1][1]*self.w[1][2]*self.w[3][1]/1080 + x**2*self.w[1][1]*self.w[2][1]**2/270 + x**2*self.w[1][1]*self.w[4][2]/72 + x**2*self.w[2][1]*self.w[3][1]/24 - x**2*self.w[5][3]/36 - x*self.w[1][1]**3*self.w[1][2]**3/90720 - x*self.w[1][1]**2*self.w[1][2]*self.w[3][2]/360 + 7*x*self.w[1][1]*self.w[1][2]**2*self.w[3][1]/2160 - x*self.w[1][1]*self.w[1][2]*self.w[2][1]**2/540 + x*self.w[1][1]*self.w[5][5]/12 + x*self.w[1][2]*self.w[5][3]/36 - x*self.w[2][1]*self.w[4][2]/36 + x*self.w[3][1]*self.w[3][2]/9 + x*self.w[6][6]/2 + self.w[7][12]
                bet[7][13] = -x**4*self.w[1][1]**3/1260 - 17*x**3*self.w[1][1]**3*self.w[1][2]/7560 - x**3*self.w[1][1]*self.w[3][1]/90 - 13*x**2*self.w[1][1]**3*self.w[1][2]**2/10080 - x**2*self.w[1][1]**2*self.w[3][2]/60 - x**2*self.w[1][1]*self.w[1][2]*self.w[3][1]/30 - x**2*self.w[1][1]*self.w[2][1]**2/120 - x*self.w[1][1]**3*self.w[1][2]**3/6048 - x*self.w[1][1]**2*self.w[1][2]*self.w[3][2]/180 - x*self.w[1][1]*self.w[1][2]**2*self.w[3][1]/180 - x*self.w[1][1]*self.w[1][2]*self.w[2][1]**2/360 + x*self.w[1][1]*self.w[5][4]/3 - x*self.w[1][1]*self.w[5][5]/6 + x*self.w[2][1]*self.w[4][2]/3 - x*self.w[3][1]*self.w[3][2]/3 + self.w[7][13]
                bet[7][14] = x**5*self.w[1][1]**2/5040 + 23*x**4*self.w[1][1]**2*self.w[1][2]/30240 + x**4*self.w[3][1]/360 - 17*x**3*self.w[1][1]**2*self.w[1][2]**2/30240 - x**3*self.w[1][1]*self.w[1][2]*self.w[2][1]/144 - x**3*self.w[1][1]*self.w[3][2]/360 + x**3*self.w[1][2]*self.w[3][1]/90 - x**3*self.w[2][1]**2/120 - x**2*self.w[1][1]**2*self.w[1][2]**3/1260 - x**2*self.w[1][1]*self.w[1][2]**2*self.w[2][1]/288 - x**2*self.w[1][1]*self.w[1][2]*self.w[3][2]/120 + x**2*self.w[1][1]*self.w[4][3]/24 - x**2*self.w[1][2]**2*self.w[3][1]/90 + x**2*self.w[1][2]*self.w[2][1]**2/360 + x**2*self.w[1][2]*self.w[4][2]/12 - x**2*self.w[2][1]*self.w[3][2]/24 - x**2*self.w[5][4]/6 + x**2*self.w[5][5]/12 - x*self.w[1][1]**2*self.w[1][2]**4/6048 - x*self.w[1][1]*self.w[1][2]**2*self.w[3][2]/720 - x*self.w[1][2]**3*self.w[3][1]/360 + x*self.w[1][2]**2*self.w[2][1]**2/360 + x*self.w[1][2]*self.w[5][4]/6 - x*self.w[1][2]*self.w[5][5]/12 - x*self.w[2][1]*self.w[4][3]/12 - x*self.w[3][2]**2/12 + x*self.w[6][7]/2 + self.w[7][14]
                bet[7][15] = x**4*self.w[1][1]**3/3360 + 11*x**3*self.w[1][1]**3*self.w[1][2]/15120 + x**3*self.w[1][1]*self.w[3][1]/180 + x**2*self.w[1][1]**3*self.w[1][2]**2/2520 + x**2*self.w[1][1]*self.w[1][2]*self.w[3][1]/60 + x*self.w[1][1]**3*self.w[1][2]**3/30240 + x*self.w[1][1]*self.w[1][2]**2*self.w[3][1]/360 - x*self.w[1][1]*self.w[5][4]/6 - x*self.w[2][1]*self.w[4][2]/6 + x*self.w[3][1]*self.w[3][2]/6 + self.w[7][15]
                bet[7][16] = x**5*self.w[1][1]**2/10080 + x**4*self.w[1][1]**2*self.w[1][2]/3780 - x**4*self.w[1][1]*self.w[2][1]/1440 - x**4*self.w[3][1]/720 + x**3*self.w[1][1]**2*self.w[1][2]**2/1890 + x**3*self.w[1][1]*self.w[1][2]*self.w[2][1]/720 + x**3*self.w[1][1]*self.w[3][2]/72 - x**3*self.w[1][2]*self.w[3][1]/180 + x**3*self.w[2][1]**2/90 + x**2*self.w[1][1]**2*self.w[1][2]**3/3360 + x**2*self.w[1][1]*self.w[1][2]**2*self.w[2][1]/720 + x**2*self.w[1][1]*self.w[1][2]*self.w[3][2]/72 - x**2*self.w[1][1]*self.w[4][3]/12 + x**2*self.w[1][2]**2*self.w[3][1]/180 + x**2*self.w[1][2]*self.w[2][1]**2/180 - x**2*self.w[1][2]*self.w[4][2]/24 + x**2*self.w[2][1]*self.w[3][2]/24 + x**2*self.w[5][4]/12 + x*self.w[1][1]**2*self.w[1][2]**4/30240 + x*self.w[1][1]*self.w[5][6]/12 + x*self.w[1][2]**3*self.w[3][1]/720 - x*self.w[1][2]**2*self.w[2][1]**2/720 - x*self.w[1][2]*self.w[5][4]/12 + x*self.w[2][1]*self.w[4][3]/4 + x*self.w[3][2]**2/6 + x*self.w[6][8]/2 + self.w[7][16]
                bet[7][17] = -x**5*self.w[1][1]**2/10080 - x**4*self.w[1][1]**2*self.w[1][2]/2016 - 23*x**3*self.w[1][1]**2*self.w[1][2]**2/30240 - x**3*self.w[1][1]*self.w[3][2]/180 - x**3*self.w[2][1]**2/360 - x**2*self.w[1][1]**2*self.w[1][2]**3/3360 - x**2*self.w[1][1]*self.w[1][2]*self.w[3][2]/60 - x**2*self.w[1][2]*self.w[2][1]**2/120 - x*self.w[1][1]**2*self.w[1][2]**4/30240 - x*self.w[1][1]*self.w[1][2]**2*self.w[3][2]/360 - x*self.w[1][1]*self.w[5][6]/6 - x*self.w[1][2]**2*self.w[2][1]**2/720 - x*self.w[2][1]*self.w[4][3]/6 - x*self.w[3][2]**2/12 + self.w[7][17]
                bet[7][18] = x**6*self.w[1][1]/30240 + x**5*self.w[1][1]*self.w[1][2]/5040 + x**4*self.w[1][1]*self.w[1][2]**2/3780 - x**4*self.w[1][2]*self.w[2][1]/1440 + x**4*self.w[3][2]/720 - x**3*self.w[1][1]*self.w[1][2]**3/3780 - x**3*self.w[1][2]**2*self.w[2][1]/360 + x**3*self.w[1][2]*self.w[3][2]/180 - x**2*self.w[1][1]*self.w[1][2]**4/5040 - x**2*self.w[1][2]**3*self.w[2][1]/1440 - x**2*self.w[1][2]**2*self.w[3][2]/180 - x**2*self.w[1][2]*self.w[4][3]/24 + x**2*self.w[5][6]/12 - x*self.w[1][1]*self.w[1][2]**5/30240 - x*self.w[1][2]**3*self.w[3][2]/720 - x*self.w[1][2]*self.w[5][6]/12 + x*self.w[6][9]/2 + self.w[7][18]                
        else:
            bet[1][1] = x + self.w[1][1]
            bet[1][2] = self.w[1][2]
            bet[2][1] = x*self.w[1][2]/2 + self.w[2][1]
            bet[3][1] = x**2*self.w[1][2]/12 - x*self.w[1][1]*self.w[1][2]/12 + x*self.w[2][1]/2 + self.w[3][1]
            bet[3][2] = -x*self.w[1][2]**2/12 + self.w[3][2]
            bet[4][1] = -x**2*self.w[1][1]*self.w[1][2]/24 + x**2*self.w[2][1]/12 - x*self.w[1][1]*self.w[2][1]/12 + x*self.w[3][1]/2 + self.w[4][1]
            bet[4][2] = -x**2*self.w[1][2]**2/24 - x*self.w[1][2]*self.w[2][1]/12 + x*self.w[3][2]/2 + self.w[4][2]
            bet[4][3] = self.w[4][3]
            if self.depth >= 5:
                bet[5][1] = -x**4*self.w[1][2]/720 - x**3*self.w[1][1]*self.w[1][2]/180 + x**2*self.w[1][1]**2*self.w[1][2]/180 - x**2*self.w[1][1]*self.w[2][1]/24 + x**2*self.w[3][1]/12 + x*self.w[1][1]**3*self.w[1][2]/720 - x*self.w[1][1]*self.w[3][1]/12 + x*self.w[4][1]/2 + self.w[5][1]
                bet[5][2] = x**3*self.w[1][2]**2/360 + x**2*self.w[1][1]*self.w[1][2]**2/120 + x*self.w[1][1]**2*self.w[1][2]**2/720 - x*self.w[1][2]*self.w[3][1]/6 + x*self.w[2][1]**2/12 + self.w[5][2]
                bet[5][3] = x**3*self.w[1][2]**2/120 - x**2*self.w[1][1]*self.w[1][2]**2/360 + x**2*self.w[1][2]*self.w[2][1]/24 - x**2*self.w[3][2]/12 - x*self.w[1][1]**2*self.w[1][2]**2/360 + x*self.w[1][1]*self.w[3][2]/12 - x*self.w[1][2]*self.w[3][1]/12 + x*self.w[2][1]**2/12 - x*self.w[4][2]/2 + self.w[5][3]
                bet[5][4] = x**2*self.w[1][2]**3/120 + x*self.w[1][1]*self.w[1][2]**3/360 - x*self.w[1][2]*self.w[3][2]/6 + self.w[5][4]
                bet[5][5] = x**2*self.w[1][2]**3/360 - x*self.w[1][1]*self.w[1][2]**3/720 - x*self.w[1][2]*self.w[3][2]/12 + x*self.w[4][3]/2 + self.w[5][5]
                bet[5][6] = -x*self.w[1][2]**4/720 + self.w[5][6]
            if self.depth >= 6:
                bet[6][1] = x**4*self.w[1][1]*self.w[1][2]/1440 - x**4*self.w[2][1]/720 + x**3*self.w[1][1]**2*self.w[1][2]/360 - x**3*self.w[1][1]*self.w[2][1]/180 + x**2*self.w[1][1]**3*self.w[1][2]/1440 + x**2*self.w[1][1]**2*self.w[2][1]/180 - x**2*self.w[1][1]*self.w[3][1]/24 + x**2*self.w[4][1]/12 + x*self.w[1][1]**3*self.w[2][1]/720 - x*self.w[1][1]*self.w[4][1]/12 + x*self.w[5][1]/2 + self.w[6][1]
                bet[6][2] = x**4*self.w[1][2]**2/1440 - x**3*self.w[1][1]*self.w[1][2]**2/720 + x**3*self.w[1][2]*self.w[2][1]/120 - x**2*self.w[1][1]**2*self.w[1][2]**2/720 - x**2*self.w[1][1]*self.w[1][2]*self.w[2][1]/360 + x**2*self.w[1][1]*self.w[3][2]/24 - x**2*self.w[1][2]*self.w[3][1]/24 + x**2*self.w[2][1]**2/24 - x**2*self.w[4][2]/12 - x*self.w[1][1]**2*self.w[1][2]*self.w[2][1]/360 + x*self.w[1][1]*self.w[4][2]/12 - x*self.w[1][2]*self.w[4][1]/6 + x*self.w[2][1]*self.w[3][1]/12 + x*self.w[5][3]/2 + self.w[6][2]
                bet[6][3] = x**3*self.w[1][1]*self.w[1][2]**2/144 - x**3*self.w[1][2]*self.w[2][1]/72 + x**2*self.w[1][1]**2*self.w[1][2]**2/288 + x**2*self.w[1][1]*self.w[1][2]*self.w[2][1]/72 - x**2*self.w[1][1]*self.w[3][2]/12 - x**2*self.w[2][1]**2/24 + x**2*self.w[4][2]/6 + x*self.w[1][1]**2*self.w[1][2]*self.w[2][1]/144 - x*self.w[1][1]*self.w[4][2]/6 + x*self.w[1][2]*self.w[4][1]/12 - x*self.w[2][1]*self.w[3][1]/12 + x*self.w[5][2]/2 - x*self.w[5][3] + self.w[6][3]
                bet[6][4] = x**3*self.w[1][2]**3/240 + x**2*self.w[1][1]*self.w[1][2]**3/720 + x**2*self.w[1][2]**2*self.w[2][1]/120 - x**2*self.w[1][2]*self.w[3][2]/12 + x*self.w[1][1]*self.w[1][2]**2*self.w[2][1]/360 - x*self.w[1][2]*self.w[4][2]/12 - x*self.w[2][1]*self.w[3][2]/12 + x*self.w[5][4]/2 + self.w[6][4]
                bet[6][5] = -x*self.w[1][2]*self.w[4][2]/18 + x*self.w[2][1]*self.w[3][2]/18 + self.w[6][5]
                bet[6][6] = x**3*self.w[1][2]**3/720 - x**2*self.w[1][1]*self.w[1][2]**3/1440 + x**2*self.w[1][2]**2*self.w[2][1]/360 - x**2*self.w[1][2]*self.w[3][2]/24 + x**2*self.w[4][3]/12 - x*self.w[1][1]*self.w[1][2]**2*self.w[2][1]/720 - x*self.w[1][1]*self.w[4][3]/12 - x*self.w[1][2]*self.w[4][2]/18 - x*self.w[2][1]*self.w[3][2]/36 + x*self.w[5][5]/2 + self.w[6][6]
                bet[6][7] = -x*self.w[1][2]*self.w[4][3]/6 + self.w[6][7]
                bet[6][8] = -x**2*self.w[1][2]**4/1440 - x*self.w[1][2]**3*self.w[2][1]/720 + x*self.w[1][2]*self.w[4][3]/12 + x*self.w[5][6]/2 + self.w[6][8]
                bet[6][9] = self.w[6][9]
            if self.depth >= 7:
                bet[7][1] = x**6*self.w[1][2]/30240 + x**5*self.w[1][1]*self.w[1][2]/5040 + x**4*self.w[1][1]**2*self.w[1][2]/3780 + x**4*self.w[1][1]*self.w[2][1]/1440 - x**4*self.w[3][1]/720 - x**3*self.w[1][1]**3*self.w[1][2]/3780 + x**3*self.w[1][1]**2*self.w[2][1]/360 - x**3*self.w[1][1]*self.w[3][1]/180 - x**2*self.w[1][1]**4*self.w[1][2]/5040 + x**2*self.w[1][1]**3*self.w[2][1]/1440 + x**2*self.w[1][1]**2*self.w[3][1]/180 - x**2*self.w[1][1]*self.w[4][1]/24 + x**2*self.w[5][1]/12 - x*self.w[1][1]**5*self.w[1][2]/30240 + x*self.w[1][1]**3*self.w[3][1]/720 - x*self.w[1][1]*self.w[5][1]/12 + x*self.w[6][1]/2 + self.w[7][1]
                bet[7][2] = -x**5*self.w[1][2]**2/10080 - x**4*self.w[1][1]*self.w[1][2]**2/2016 - 23*x**3*self.w[1][1]**2*self.w[1][2]**2/30240 + x**3*self.w[1][2]*self.w[3][1]/180 - x**3*self.w[2][1]**2/360 - x**2*self.w[1][1]**3*self.w[1][2]**2/3360 + x**2*self.w[1][1]*self.w[1][2]*self.w[3][1]/60 - x**2*self.w[1][1]*self.w[2][1]**2/120 - x*self.w[1][1]**4*self.w[1][2]**2/30240 + x*self.w[1][1]**2*self.w[1][2]*self.w[3][1]/360 - x*self.w[1][1]**2*self.w[2][1]**2/720 - x*self.w[1][2]*self.w[5][1]/6 + x*self.w[2][1]*self.w[4][1]/6 - x*self.w[3][1]**2/12 + self.w[7][2]
                bet[7][3] = x**5*self.w[1][2]**2/10080 + x**4*self.w[1][1]*self.w[1][2]**2/3780 + x**4*self.w[1][2]*self.w[2][1]/1440 + x**4*self.w[3][2]/720 + x**3*self.w[1][1]**2*self.w[1][2]**2/1890 - x**3*self.w[1][1]*self.w[1][2]*self.w[2][1]/720 + x**3*self.w[1][1]*self.w[3][2]/180 - x**3*self.w[1][2]*self.w[3][1]/72 + x**3*self.w[2][1]**2/90 + x**2*self.w[1][1]**3*self.w[1][2]**2/3360 - x**2*self.w[1][1]**2*self.w[1][2]*self.w[2][1]/720 - x**2*self.w[1][1]**2*self.w[3][2]/180 - x**2*self.w[1][1]*self.w[1][2]*self.w[3][1]/72 + x**2*self.w[1][1]*self.w[2][1]**2/180 + x**2*self.w[1][1]*self.w[4][2]/24 - x**2*self.w[1][2]*self.w[4][1]/12 + x**2*self.w[2][1]*self.w[3][1]/24 + x**2*self.w[5][3]/12 + x*self.w[1][1]**4*self.w[1][2]**2/30240 - x*self.w[1][1]**3*self.w[3][2]/720 - x*self.w[1][1]**2*self.w[2][1]**2/720 - x*self.w[1][1]*self.w[5][3]/12 + x*self.w[1][2]*self.w[5][1]/12 - x*self.w[2][1]*self.w[4][1]/4 + x*self.w[3][1]**2/6 + x*self.w[6][2]/2 + self.w[7][3]
                bet[7][4] = x**4*self.w[1][2]**3/3360 + 11*x**3*self.w[1][1]*self.w[1][2]**3/15120 - x**3*self.w[1][2]*self.w[3][2]/180 + x**2*self.w[1][1]**2*self.w[1][2]**3/2520 - x**2*self.w[1][1]*self.w[1][2]*self.w[3][2]/60 + x*self.w[1][1]**3*self.w[1][2]**3/30240 - x*self.w[1][1]**2*self.w[1][2]*self.w[3][2]/360 - x*self.w[1][2]*self.w[5][3]/6 - x*self.w[2][1]*self.w[4][2]/6 + x*self.w[3][1]*self.w[3][2]/6 + self.w[7][4]
                bet[7][5] = x**5*self.w[1][2]**2/5040 + 23*x**4*self.w[1][1]*self.w[1][2]**2/30240 - x**4*self.w[3][2]/360 - 17*x**3*self.w[1][1]**2*self.w[1][2]**2/30240 + x**3*self.w[1][1]*self.w[1][2]*self.w[2][1]/144 - x**3*self.w[1][1]*self.w[3][2]/90 + x**3*self.w[1][2]*self.w[3][1]/360 - x**3*self.w[2][1]**2/120 - x**2*self.w[1][1]**3*self.w[1][2]**2/1260 + x**2*self.w[1][1]**2*self.w[1][2]*self.w[2][1]/288 + x**2*self.w[1][1]**2*self.w[3][2]/90 + x**2*self.w[1][1]*self.w[1][2]*self.w[3][1]/120 + x**2*self.w[1][1]*self.w[2][1]**2/360 - x**2*self.w[1][1]*self.w[4][2]/12 + x**2*self.w[1][2]*self.w[4][1]/24 - x**2*self.w[2][1]*self.w[3][1]/24 + x**2*self.w[5][2]/12 - x**2*self.w[5][3]/6 - x*self.w[1][1]**4*self.w[1][2]**2/6048 + x*self.w[1][1]**3*self.w[3][2]/360 + x*self.w[1][1]**2*self.w[1][2]*self.w[3][1]/720 + x*self.w[1][1]**2*self.w[2][1]**2/360 - x*self.w[1][1]*self.w[5][2]/12 + x*self.w[1][1]*self.w[5][3]/6 + x*self.w[2][1]*self.w[4][1]/12 - x*self.w[3][1]**2/12 + x*self.w[6][3]/2 + self.w[7][5]
                bet[7][6] = -x**4*self.w[1][2]**3/1260 - 17*x**3*self.w[1][1]*self.w[1][2]**3/7560 + x**3*self.w[1][2]*self.w[3][2]/90 - 13*x**2*self.w[1][1]**2*self.w[1][2]**3/10080 + x**2*self.w[1][1]*self.w[1][2]*self.w[3][2]/30 + x**2*self.w[1][2]**2*self.w[3][1]/60 - x**2*self.w[1][2]*self.w[2][1]**2/120 - x*self.w[1][1]**3*self.w[1][2]**3/6048 + x*self.w[1][1]**2*self.w[1][2]*self.w[3][2]/180 + x*self.w[1][1]*self.w[1][2]**2*self.w[3][1]/180 - x*self.w[1][1]*self.w[1][2]*self.w[2][1]**2/360 - x*self.w[1][2]*self.w[5][2]/6 + x*self.w[1][2]*self.w[5][3]/3 + x*self.w[2][1]*self.w[4][2]/3 - x*self.w[3][1]*self.w[3][2]/3 + self.w[7][6]
                bet[7][7] = x**4*self.w[1][2]**3/1120 + x**3*self.w[1][1]*self.w[1][2]**3/10080 + x**3*self.w[1][2]**2*self.w[2][1]/240 - x**3*self.w[1][2]*self.w[3][2]/60 - x**2*self.w[1][1]**2*self.w[1][2]**3/5040 + x**2*self.w[1][1]*self.w[1][2]**2*self.w[2][1]/720 + x**2*self.w[1][1]*self.w[1][2]*self.w[3][2]/180 - x**2*self.w[1][2]**2*self.w[3][1]/360 + x**2*self.w[1][2]*self.w[2][1]**2/180 - x**2*self.w[1][2]*self.w[4][2]/24 - x**2*self.w[2][1]*self.w[3][2]/24 + x**2*self.w[5][4]/12 - x*self.w[1][1]**3*self.w[1][2]**3/7560 + x*self.w[1][1]**2*self.w[1][2]*self.w[3][2]/180 - x*self.w[1][1]*self.w[1][2]**2*self.w[3][1]/180 + x*self.w[1][1]*self.w[1][2]*self.w[2][1]**2/240 - x*self.w[1][1]*self.w[5][4]/12 - x*self.w[1][2]*self.w[5][3]/12 - x*self.w[2][1]*self.w[4][2]/6 + x*self.w[3][1]*self.w[3][2]/12 + x*self.w[6][4]/2 + self.w[7][7]
                bet[7][8] = -x**3*self.w[1][2]**4/1120 - x**2*self.w[1][1]*self.w[1][2]**4/1680 + x**2*self.w[1][2]**2*self.w[3][2]/40 - x*self.w[1][1]**2*self.w[1][2]**4/7560 + x*self.w[1][1]*self.w[1][2]**2*self.w[3][2]/120 - x*self.w[1][2]**3*self.w[3][1]/360 + x*self.w[1][2]**2*self.w[2][1]**2/720 - x*self.w[1][2]*self.w[5][4]/6 - x*self.w[3][2]**2/12 + self.w[7][8]
                bet[7][9] = x**4*self.w[1][2]**3/5040 + 17*x**3*self.w[1][1]*self.w[1][2]**3/30240 - x**3*self.w[1][2]*self.w[3][2]/360 + x**2*self.w[1][1]**2*self.w[1][2]**3/30240 - x**2*self.w[1][1]*self.w[1][2]*self.w[3][2]/120 - x**2*self.w[1][2]**2*self.w[3][1]/90 + x**2*self.w[1][2]*self.w[2][1]**2/180 - x**2*self.w[1][2]*self.w[4][2]/36 + x**2*self.w[2][1]*self.w[3][2]/36 - x*self.w[1][1]**3*self.w[1][2]**3/18144 - x*self.w[1][1]**2*self.w[1][2]*self.w[3][2]/720 + x*self.w[1][1]*self.w[1][2]**2*self.w[3][1]/1080 - x*self.w[1][1]*self.w[1][2]*self.w[2][1]**2/2160 + x*self.w[1][2]*self.w[5][2]/12 - x*self.w[1][2]*self.w[5][3]/36 - x*self.w[2][1]*self.w[4][2]/12 + 5*x*self.w[3][1]*self.w[3][2]/36 + x*self.w[6][5]/2 + self.w[7][9]
                bet[7][10] = x**3*self.w[1][2]**4/30240 - x**2*self.w[1][1]*self.w[1][2]**4/7560 - x**2*self.w[1][2]**2*self.w[3][2]/360 - x*self.w[1][1]**2*self.w[1][2]**4/18144 - x*self.w[1][1]*self.w[1][2]**2*self.w[3][2]/1080 + x*self.w[1][2]**3*self.w[3][1]/360 - x*self.w[1][2]**2*self.w[2][1]**2/720 + x*self.w[3][2]**2/36 + self.w[7][10]
                bet[7][11] = x**4*self.w[1][2]**3/3024 - x**3*self.w[1][1]*self.w[1][2]**3/15120 + x**3*self.w[1][2]**2*self.w[2][1]/720 - x**3*self.w[1][2]*self.w[3][2]/120 + x**2*self.w[1][1]**2*self.w[1][2]**3/7560 - x**2*self.w[1][1]*self.w[1][2]**2*self.w[2][1]/1440 + x**2*self.w[1][1]*self.w[1][2]*self.w[3][2]/360 - x**2*self.w[1][1]*self.w[4][3]/24 - x**2*self.w[1][2]**2*self.w[3][1]/360 + x**2*self.w[1][2]*self.w[2][1]**2/360 - x**2*self.w[1][2]*self.w[4][2]/36 - x**2*self.w[2][1]*self.w[3][2]/72 + x**2*self.w[5][5]/12 + x*self.w[1][1]**3*self.w[1][2]**3/90720 + x*self.w[1][1]**2*self.w[1][2]*self.w[3][2]/360 - 7*x*self.w[1][1]*self.w[1][2]**2*self.w[3][1]/2160 + x*self.w[1][1]*self.w[1][2]*self.w[2][1]**2/1080 - x*self.w[1][1]*self.w[5][5]/12 - x*self.w[1][2]*self.w[5][3]/36 - x*self.w[2][1]*self.w[4][2]/12 + x*self.w[3][1]*self.w[3][2]/18 + x*self.w[6][6]/2 + self.w[7][11]
                bet[7][12] = -x**3*self.w[1][2]**4/10080 + x**2*self.w[1][1]*self.w[1][2]**4/6048 + x**2*self.w[1][2]**2*self.w[3][2]/180 + x*self.w[1][1]**2*self.w[1][2]**4/90720 + x*self.w[1][1]*self.w[1][2]**2*self.w[3][2]/540 - x*self.w[1][2]**3*self.w[3][1]/360 + x*self.w[1][2]**2*self.w[2][1]**2/720 - x*self.w[1][2]*self.w[5][5]/6 + x*self.w[2][1]*self.w[4][3]/6 + x*self.w[3][2]**2/36 + self.w[7][12]
                bet[7][13] = -x**3*self.w[1][2]**4/1260 - x**2*self.w[1][1]*self.w[1][2]**4/3360 + x**2*self.w[1][2]**2*self.w[3][2]/40 - x**2*self.w[1][2]*self.w[4][3]/12 + x*self.w[1][1]**2*self.w[1][2]**4/6048 + x*self.w[1][1]*self.w[1][2]**2*self.w[3][2]/720 + x*self.w[1][2]**3*self.w[3][1]/360 - x*self.w[1][2]**2*self.w[2][1]**2/720 - x*self.w[1][2]*self.w[5][4]/6 + x*self.w[1][2]*self.w[5][5]/12 - x*self.w[2][1]*self.w[4][3]/4 - x*self.w[3][2]**2/6 + x*self.w[6][7]/2 + self.w[7][13]
                bet[7][14] = x**2*self.w[1][2]**5/5040 + x*self.w[1][1]*self.w[1][2]**5/6048 - x*self.w[1][2]**3*self.w[3][2]/360 + self.w[7][14]
                bet[7][15] = x**3*self.w[1][2]**4/3360 + x**2*self.w[1][1]*self.w[1][2]**4/5040 - x**2*self.w[1][2]**3*self.w[2][1]/1440 - x**2*self.w[1][2]**2*self.w[3][2]/72 + x**2*self.w[1][2]*self.w[4][3]/24 + x**2*self.w[5][6]/12 - x*self.w[1][1]**2*self.w[1][2]**4/30240 - x*self.w[1][1]*self.w[5][6]/12 - x*self.w[1][2]**3*self.w[3][1]/720 + x*self.w[1][2]*self.w[5][4]/12 + x*self.w[2][1]*self.w[4][3]/12 + x*self.w[3][2]**2/12 + x*self.w[6][8]/2 + self.w[7][15]
                bet[7][16] = x**2*self.w[1][2]**5/10080 - x*self.w[1][1]*self.w[1][2]**5/30240 - x*self.w[1][2]*self.w[5][6]/6 + self.w[7][16]
                bet[7][17] = -x**2*self.w[1][2]**5/10080 + x*self.w[1][1]*self.w[1][2]**5/30240 + x*self.w[1][2]**3*self.w[3][2]/720 + x*self.w[1][2]*self.w[5][6]/12 + x*self.w[6][9]/2 + self.w[7][17]
                bet[7][18] = x*self.w[1][2]**6/30240 + self.w[7][18]                
        for i in range(len(bet)):
            for j in range(len(bet[i])):
                bet[i][j] = bet[i][j].expand()
        self.w = bet.copy()

def printd(text):
    print(colored("DEBUG:", "yellow", attrs=['bold']), text)
    
def printe(text):
    print(colored("ERROR:", "red", attrs=['bold']), text)

def bch7(A, B, depth = 7):
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
        e51 = Corxet(A, e41)
        e52 = Corxet(B, e41)
        e53 = Corxet(A, -e42)
        e54 = Corxet(B, e42)
        e55 = Corxet(A, e43)
        e56 = Corxet(B, e43)
        D += sp.Rational(1, 720)*(-e51 - e56 + sp.S(6)*e53 + sp.S(6)*e54 + sp.S(2)*e55 + sp.S(2)*e52)
    if (depth >= 6):
        e62 = Corxet(B, e51)
        e64 = Corxet(A, e54)
        e66 = Corxet(A, e55)
        e68 = Corxet(A, e56)
        D += sp.Rational(1, 1440)*(e62 - e68 + sp.S(2)*e66 + sp.S(6)*e64)
    if (depth >= 7):
        cc = sp.S([1, 6, -3, 3, 30, 6, -4, 9, 3, -3, -24, -18, 9, -9, -12, 4, 3, -1])
        c = []
        ca1 = Corxet(B, Corxet(A, e42))
        ca2 = Corxet(B, e54)
        e61 = Corxet(A, e51)
        e63 = Corxet(A, e52)
        e65 = Corxet(B, e52)
        e67 = Corxet(B, e55)
        e69 = Corxet(B, e56)
        c.append(Corxet(A, e61))
        c.append(Corxet(A, e63))
        c.append(Corxet(A, e64))
        c.append(Corxet(A, e62))
        c.append(Corxet(A, ca1))
        c.append(Corxet(A, -e67))
        c.append(Corxet(A, e65))
        c.append(Corxet(A, ca2))
        c.append(Corxet(A, -e69))
        c.append(Corxet(B, e61))
        c.append(Corxet(B, e63))
        c.append(Corxet(B, e64))
        c.append(Corxet(B, e62))
        c.append(Corxet(B, ca1))
        c.append(Corxet(B, -e67))
        c.append(Corxet(B, e65))
        c.append(Corxet(B, ca2))
        c.append(Corxet(B, -e69))
        c7 = sp.S(0)
        for i in range(0, len(c)):
            c7 += cc[i]*c[i]
        D += sp.Rational(1, 30240)*c7
    return D
