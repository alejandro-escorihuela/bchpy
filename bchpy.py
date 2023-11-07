#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 01-12-2020
# alex
# bchpy.py

import numpy as np
import sympy as sp
import relations as rl
import recursive as rc
import bch20Lyndon as b20
import time as tm
import datetime
import sys
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
 
    def __init__(self, i = 1, j = 1, basis_type = "!"):
        self.i = i
        self.j = j
        if basis_type not in ["E", "C", "M", "M4", "M6"]:
            printe("%s no és una base vàlida" % (basis_type))
            exit(-1)
        self.t = basis_type
        if i == 0 or j == 0:
            self.is_zero = True
    
    def _sympystr(self, printer, *args):
        if self.is_zero:
            return "0"
        return "%s%s%s" % (self.t, printer._print(self.args[0]), printer._print(self.args[1]))
    
    def _sympyrepr(self, printer, *args):
        nom = self.__class__.__name__
        return "%s(%d,%d,'%s')" % (nom, self.i, self.j, self.t)
    
    def _pretty(self, printer, *args):
        if self.is_zero:
            return prettyForm(pretty_symbol("0"))
        return prettyForm(pretty_symbol("%s_%s%s" % (self.t, printer._print(self.args[0]), printer._print(self.args[1]))))
    def _latex(self, printer, *args):
        return "%s_{%d,%d}" % (self.t, self.i, self.j)         
    
class Claudator(sp.Expr):
    def __new__(cls, A, B):
        return cls.eval(cls, A, B)

    def eval(cls, A, B):
        if not A != B:
            return sp.S.Zero
        if A.is_commutative or B.is_commutative:
            return sp.S.Zero
        if isinstance(A, sp.Add):
            sargs = []
            for term in A.args:
                sargs.append(Claudator(term, B))
            return sp.Add(*sargs)
        if isinstance(B, sp.Add):
            sargs = []
            for term in B.args:
                sargs.append(Claudator(A, term))
            return sp.Add(*sargs)
        cA, ncA = A.args_cnc()
        cB, ncB = B.args_cnc()
        c_part = cA + cB
        if c_part:
            return sp.Mul(sp.Mul(*c_part), cls(sp.Mul._from_args(ncA), sp.Mul._from_args(ncB)))
        if A.compare(B) == 1:
            return sp.S.NegativeOne*cls.eval(cls, B, A)
        i, j, k, l = A.i, A.j, B.i, B.j
        if A.t != B.t:
            printe("Els objectes %s i %s no pertanyen a la mateixa base." % (A, B))
            exit(-1)
        if A.t == "E":
            ll = rl.corxetErel(i, j, k, l)
        elif A.t == "C":
            ll = rl.corxetCrel(i, j, k, l)
        elif A.t == "M":
            ll = rl.corxetMrel(i, j, k, l)
        elif A.t == "M4":
            ll = rl.corxetM4rel(i, j, k, l)
        elif A.t == "M6":
            ll = rl.corxetM6rel(i, j, k, l)            
        else:
            printe("No hi ha cap base definida amb %s." % (A.t))
            exit(-1)            
        args = []
        for it in ll:
            elem = sp.Mul(sp.Rational(it[0], it[1]), Eel(it[2], it[3], A.t))
            args.append(elem)
        return sp.Add(*args)

class Metode():
    depth = 0
    w = []
    setw = False
    num = False
    rknd = False # Only for symmetrical iteration
    cofs = []
    t = ""
    
    def __init__(self, depth = 6, basis_type = "E", numeric = False, rkn = False):
        self.num = numeric
        self.rknd = rkn
        dict_depth = {
            "E"  : [4, len(rl.tamE)],
            "C"  : [4, len(rl.tamC)],
            "M"  : [5, len(rl.tamM)],
            "M4" : [7, len(rl.tamM4)],
            "M6" : [9, len(rl.tamM6)]          
        }
        for it in dict_depth:
            if basis_type == it:
                self.t = it
        if self.t == "":
            printe("%s no és una base vàlida" % (basis_type))
            exit(-1)            
        if depth < dict_depth[self.t][0]:
            self.depth = dict_depth[self.t][0]
        elif depth >= dict_depth[self.t][1]:
            self.depth = dict_depth[self.t][1]
        else:
            self.depth = depth
        self.w = self.__init_mat()
        
    def setBAB(self, *args, debug = False):
        if self.t != "E":
            printe("Per mètodes BAB la base ha de ser E i no %s" % (self.t))
            exit(-1)            
        if self.setw == True:
            self.w = self.__init_mat()
        self.cofs = list(args)
        # cof_li = list(reversed(args))
        cof_li = self.cofs
        tipus = palindromic(cof_li)
        senar = (len(cof_li) % 2) != 0
        if tipus == 1 and senar == True:
            cof_li = cof_li[len(cof_li)//2:]
            iniEl = 0 if len(cof_li)%2 == 0 else 1
            self.w[1][iniEl + 1] = cof_li[0]
            for i in range(1, len(cof_li)):
                if debug == True:
                    txt_p = "Iteració simètrica " + str(i) + " amb Eel(1, " + str(iniEl + 1 - (i%2)) + ")"
                    printd(txt_p)
                self.__recur_AB(cof_li[i], 2 + (iniEl + i)%2)            
        else:        
            self.w[1][2] = cof_li[0]
            for i in range(1, len(cof_li)):
                if debug == True:
                    txt_p = "Iteració " + str(i) + " amb Eel(1, " + str(2 - (i%2)) + ")"
                    printd(txt_p)
                self.__recur_AB(cof_li[i], (i - 1)%2)
        self.setw = True
        
    def setABA(self, *args, debug = False):
        if self.t != "E":
            printe("Per mètodes ABA la base ha de ser E i no %s" % (self.t))
            exit(-1)         
        if self.setw == True:
            self.w = self.__init_mat()
        self.cofs = list(args)
        # cof_li = list(reversed(args))
        cof_li = self.cofs
        tipus = palindromic(cof_li)
        senar = (len(cof_li) % 2) != 0
        if tipus == 1 and senar == True:
            cof_li = cof_li[len(cof_li)//2:]
            iniEl = 1 if len(cof_li)%2 == 0 else 0
            self.w[1][iniEl + 1] = cof_li[0]
            for i in range(1, len(cof_li)):
                if debug == True:
                    txt_p = "Iteració simètrica " + str(i) + " amb Eel(1, " + str(iniEl + 1 - (i%2)) + ")"
                    printd(txt_p)
                self.__recur_AB(cof_li[i], 2 + (iniEl + i)%2)              
        else:
            self.w[1][1] = cof_li[0]
            for i in range(1, len(cof_li)):
                if debug == True:
                    printd("Iteració " + str(i) + " amb Eel(1, " + str((i%2) + 1) + ")")
                self.__recur_AB(cof_li[i], i%2)             
        self.setw = True

    def setXX(self, *args, debug = False):
        if self.t != "C":
            printe("Per mètodes XX* la base ha de ser Z i no %s" % (self.t))
            exit(-1)         
        if self.setw == True:
            self.w = self.__init_mat()        
        self.cofs = list(args)
        #cofsR = list(reversed(args))
        for i in range(1, self.depth + 1):
            self.w[i][1] = self.cofs[0]**i
        for i in range(1, len(self.cofs)):
            self.__recur_XX(self.cofs[i], i%2)
        self.setw = True

    def setSS2(self, *args, debug = False):
        return self.setSS(*args, base = 2, debug = debug)

    def setSS4(self, *args, debug = False):
        return self.setSS(*args, base = 4, debug = debug)

    def setSS6(self, *args, debug = False):
        return self.setSS(*args, base = 6, debug = debug)
        
    def setSS(self, *args, base = 2, debug = False):
        if base == 2 and self.t != "M":
            printe("Per mètodes SS(S2) la base ha de ser M i no %s" % (self.t))
            exit(-1)             
        elif base != 2 and self.t != "M" + str(base):
            printe("Per mètodes SS(S" + str(base) + ") la base ha de ser M" + str(base) + " i no %s" % (self.t))
            exit(-1)         
        if self.setw == True:
            self.w = self.__init_mat()
        self.cofs = list(args)
        cof_li = self.cofs
        tipus = palindromic(cof_li)
        senar = (len(cof_li) % 2) != 0
        if tipus == 1:
            cof_li = cof_li[len(cof_li)//2:]
        if senar == True:
            self.w[1][1] = cof_li[0]
            for i in range(base + 1, self.depth + 1, 2):
                self.w[i][1] = cof_li[0]**i
        for i in range(1, len(cof_li)):
            if debug == True:
                txt_p = "Iteració " + str(i) + " amb mètode S" + str(base) + " del tipus " + str(tipus)
                printd(txt_p)
            if base == 2:
                self.__recur_S2(cof_li[i], tipus)
            elif base == 4:
                self.__recur_S4(cof_li[i], tipus)
            elif base == 6:
                self.__recur_S6(cof_li[i], tipus)
            elif base == 8:
                self.__recur_S8(cof_li[i], tipus)
            elif base == 10:
                self.__recur_S10(cof_li[i], tipus)                    
        self.setw = True
        
    def expand(self, debug = False):        
        for i in range(1, len(self.w)):
            for j in range(1, len(self.w[i])):
                if debug == True:
                    printd("Expandint " + str(i) + ", " + str(j))
                self.w[i][j] = self.w[i][j].expand()
                
    def simplify(self, debug = False):        
        for i in range(1, len(self.w)):
            for j in range(1, len(self.w[i])):
                if debug == True:
                    printd("Simplificant " + str(i) + ", " + str(j))
                self.w[i][j] = self.w[i][j].simplify()
                
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
                esq += self.w[i][j]*Eel(i, j, self.t)
        return esq

    def importFromExpr(self, esq, debug = False):
        if self.t == "E":
            tams = rl.tamE
        elif self.t == "C":
            tams = rl.tamC
        elif self.t == "M":
            tams = rl.tamM
        elif self.t == "M4":
            tams = rl.tamM4
        elif self.t == "M6":
            tams = rl.tamM6             
        _x_diff_var = sp.Symbol("_x_diff_var")
        if debug == True:
            printd("Important esquema")
        esq_e = sp.expand(esq)
        for i in range(self.depth):
            for j in range(tams[i]):
                self.w[i + 1][j + 1] = sp.diff(esq_e.subs(Eel(i + 1, j + 1, self.t), _x_diff_var), _x_diff_var)
        self.setw = True  

    def save(self, nom_fitxer):
        f = open(nom_fitxer, "w")
        f.write(str(self.depth) + "\n")
        if not self.cofs:
            f.write(str([0]) + "\n")
        else:
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
        f.write("# " + ('%02d' % ara.day) + "-" + ('%02d' % ara.month) + "-" + str(ara.year) + "\n")
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
        
    def exportC(self, nom_fitC):
        if self.setw == False:
            printe("Les equacions del mètode no estan calculades.")
            exit(-1)        
        ara = datetime.datetime.now()
        nom_fitH = nom_fitC.replace(".c", ".h")
        with open(nom_fitH, 'w') as fH:
            fH.write("/* " + ('%02d' % ara.day) + "-" + ('%02d' % ara.month) + "-" + str(ara.year) + " */\n")
            fH.write("/* " + nom_fitH + " */\n")
            fH.write("/* Fitxer generat automàticament per bchpy. No tocar! */\n")
            fH.write("#ifndef _" + nom_fitH.replace(".", "_").upper() + "\n")
            fH.write("#define _" + nom_fitH.replace(".", "_").upper() + "\n")
            fH.write("#define POTENCIA(A, B) (pow((A), (B)))\n")
            fH.write("#include <stdio.h>\n")
            fH.write("#include <stdlib.h>\n\n")
            for i in range(1, len(self.w)):
                for j in range(1, len(self.w[i])):
                    fH.write("double w" + str(i) + str(j) + "(double *, double *, double *);\n")
            fH.write("\n#endif")
        with open(nom_fitC, 'w') as fC:
            fC.write("/* " + ('%02d' % ara.day) + "-" + ('%02d' % ara.month) + "-" + str(ara.year) + " */\n")
            fC.write("/* " + nom_fitC + " */\n")
            fC.write("/* Fitxer generat automàticament per bchpy. No tocar! */\n")
            fC.write("#include <stdio.h>\n")
            fC.write("#include <stdlib.h>\n")
            fC.write("#include \"" + nom_fitH + "\"\n\n")            
            for i in range(1, len(self.w)):
                for j in range(1, len(self.w[i])):
                    cad = str(sp.ccode(self.w[i][j]))
                    cad = strsub(r'a([0-9]*)', r'a[\1]', cad)
                    cad = strsub(r'b([0-9]*)', r'br[\1]', cad)
                    cad = strsub(r're\(br\[([0-9]*)\]\)', r'br[\1]', cad)
                    cad = strsub(r'im\(br\[([0-9]*)\]\)', r'bi[\1]', cad)
                    cad = cad.replace("I*", "")
                    cad = cad.replace("pow", "POTENCIA")
                    fC.write("double w" + str(i) + str(j) + "(double * a, double * br, double * bi) {\n")
                    fC.write("  (void) bi;\n")
                    fC.write("  return " + cad + ";\n")
                    fC.write("}\n\n")
                    
    def exportnb(self, nom_fitxer):
        if self.setw == False:
            printe("Les equacions del mètode no estan calculades.")
            exit(-1)
        ara = datetime.datetime.now()
        f = open(nom_fitxer, "w")    
        f.write("(* " + ('%02d' % ara.day) + "-" + ('%02d' % ara.month) + "-" + str(ara.year) + " *)\n")
        f.write("(* " + nom_fitxer + " *)\n")
        f.write("(* Fitxer generat automàticament per bchpy. No tocar! *)\n")
        cof_txt = str(self.cofs).replace("[", "").replace("]", "")
        f.write("(* Esquema: " + cof_txt + " *)\n\n")
        f.write("eq = {")
        for i in range(1, len(self.w)):
            for j in range(1, len(self.w[i])):
                if (self.w[i][j] != 0):
                    cad = str(self.w[i][j])
                    cad = strsub(r'a([0-9]*)', r'a\1', cad)
                    cad = strsub(r'b([0-9]*)', r'br\1', cad)
                    cad = strsub(r're\(br\[([0-9]*)\]\)', r'br\1', cad)
                    cad = strsub(r'im\(br\[([0-9]*)\]\)', r'bi\1', cad)
                    cad = cad.replace("I*", "")
                    cad = cad.replace("**", "^")
                    if i == 1 and j == 1:
                        f.write(cad + " == 1")
                    elif i == 1 and j != 1:
                        f.write(",\n\t" + cad + " == 1")
                    else:
                        f.write(",\n\t" + cad + " == 0")
        f.write("}")
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
        if self.t == "E":
            tams = rl.tamE
        elif self.t == "C":
            tams = rl.tamC
        elif self.t == "M":
            tams = rl.tamM
        elif self.t == "M4":
            tams = rl.tamM4
        elif self.t == "M6":
            tams = rl.tamM6            
        zero = sp.S(0)
        if self.num:
            zero = np.float128(0.0)
        for i in range(self.depth):
            wret.append([])
            for j in range(tams[i]):
                wret[i].append(zero)
            wret[i].insert(0, zero)
        wret.insert(0, [zero, zero])
        return wret
    
    def __str__(self):
        cad = ""
        for i in range(1, len(self.w)):
            for j in range(1, len(self.w[i])):
                cad = cad + "w(" + str(i) + "," + str(j) + ") = " + str(self.w[i][j]) + "\n"
        return cad
    
    def __recur_AB(self, x, AB):
        bet = self.__init_mat()
        if AB == 0:
            bet = rc.recA(self.w, bet, x, self.depth) 
        elif AB == 1:
            bet = rc.recB(self.w, bet, x, self.depth)
        elif AB == 2:
            bet = rc.recAsim(self.w, bet, x, self.depth, self.rknd)
        elif AB == 3:
            bet = rc.recBsim(self.w, bet, x, self.depth, self.rknd)
        # for i in range(len(bet)):
        #     for j in range(len(bet[i])):
        #         bet[i][j] = bet[i][j].expand()
        self.w = bet.copy()
        
    def __recur_XX(self, x, adj):
        bet = self.__init_mat()
        if adj == 0:
            bet = rc.recXa(self.w, bet, x, self.depth) 
        elif adj == 1:
            bet = rc.recXb(self.w, bet, x, self.depth)
        self.w = bet.copy()

    def __recur_S2(self, x, tip):
        bet = self.__init_mat()
        if tip == 1:
            bet = rc.recS2sim(self.w, bet, x, self.depth)
        else:
            bet = rc.recS2(self.w, bet, x, self.depth)
        self.w = bet.copy()

    def __recur_S4(self, x, tip):
        bet = self.__init_mat()              
        if tip == 1:
            bet = rc.recS4sim(self.w, bet, x, self.depth)
        else:
            bet = rc.recS4(self.w, bet, x, self.depth)     
        self.w = bet.copy()

    def __recur_S6(self, x, tip):
        bet = self.__init_mat()
        if tip == 0:
            printe("Per mètodes SS(S6) no estan implementades les iteracions sense cap tipus de simetria")
            exit(-1)         
        bet = rc.recS6sim(self.w, bet, x, self.depth)
        self.w = bet.copy()

    def __recur_S8(self, x, tip):
        bet = self.__init_mat()
        if tip == 0:
            printe("Per mètodes SS(S8) no estan implementades les iteracions sense cap tipus de simetria")
            exit(-1)         
        bet = rc.recS8sim(self.w, bet, x, self.depth)
        self.w = bet.copy()
        
    def __recur_S10(self, x, tip):
        bet = self.__init_mat()
        if tip == 0:
            printe("Per mètodes SS(S10) no estan implementades les iteracions sense cap tipus de simetria")
            exit(-1)         
        bet = rc.recS10sim(self.w, bet, x, self.depth)
        self.w = bet.copy()
        
def palindromic(cofs):
    # return 0
    if cofs == cofs[::-1]:
        return 1
    # cofsconj = cofs.copy()
    # for i in range(len(cofsconj)):
    #     cofsconj[i] = cofsconj[i].conjugate()
    # if cofsconj == cofs[::-1]:
    #     return 2
    return 0
        
def collectEsq(esq, basis_type = "E"):
    if basis_type == "E":
        tams = rl.tamE
    elif basis_type == "C":
        tams = rl.tamC
    elif basis_type == "M":
        tams = rl.tamM
    elif basis_type == "M4":
        tams = rl.tamM4
    elif basis_type == "M6":
        tams = rl.tamM6         
    else:
        printe("%s no és una base vàlida" % (basis_type))
        exit(-1)
    esqn = sp.S(0)
    _x_diff_var = sp.Symbol("_x_diff_var")
    esq_e = sp.expand(esq)
    for i in range(len(rl.tamE)):
        for j in range(rl.tamE[i]):
            esqn += sp.diff(esq_e.subs(Eel(i + 1, j + 1, basis_type), _x_diff_var), _x_diff_var)*Eel(i + 1, j + 1, basis_type)
    return esqn

def printi(text):
    print(colored("INFO :", "blue", attrs=['bold']), text)
        
def printd(text):
    print(colored("DEBUG:", "yellow", attrs=['bold']), text, file=sys.stderr)
    
def printe(text):
    print(colored("ERROR:", "red", attrs=['bold']), text, file=sys.stderr)

def bch20(A, B, depth = 6, debug = False):
    elem = []
    for i in range(len(b20.elemBCH)):
        elem.append(sp.S(0))
    elem[1], elem[2] = A, B
    for i in range(3, len(b20.elemBCH)):
        if b20.elemBCH[i][5] < depth or (b20.elemBCH[i][5] == depth and b20.elemBCH[i][3] != 0):
            if debug == True:
                printd("BCH: Càlcul de l'element %d" % (i))
            j, k = b20.elemBCH[i][1], b20.elemBCH[i][2]
            elem[i] = Claudator(elem[j], elem[k])
    D = A + B
    for i in range(3, len(b20.elemBCH)):
        if elem[i] != 0:
            #D += elem[i]*sp.Rational(b20.elemBCH[i][3], b20.elemBCH[i][4])
            D += elem[i]*(b20.elemBCH[i][3]/b20.elemBCH[i][4])
    return D

# deprecated
def bch9(A, B, depth = 6, debug = False):  
    if debug == True:
        printd("BCH d'ordre 1 a 4")
    e21 = Claudator(A, B).expand()
    e31 = Claudator(A, e21).expand()
    e32 = Claudator(B, e21).expand()
    e41 = Claudator(A, e31).expand()
    e42 = Claudator(A, e32).expand()
    e43 = Claudator(B, -e32).expand()
    D = A + B
    D += sp.Rational(1, 2)*e21
    D += sp.Rational(1, 12)*(e31 - e32)
    D += sp.Rational(-1, 24)*e42
    
    if (depth >= 5):
        if debug == True:
            printd("BCH d'ordre 5")
        den = sp.S([-720, -120, -360, 360, 120, 720])
        f5 = [sp.S(0)]*len(den)
        f5[0] = Claudator(A, e41)
        f5[1] = Claudator(A, e42)
        f5[2] = Claudator(A, -e43)
        f5[3] = Claudator(B, e41)
        f5[4] = Claudator(B, e42)
        f5[5] = Claudator(B, -e43)
        c5 = sp.S(0)
        for i in range(0, len(den)):
            f5[i] = f5[i].expand()
            c5 += sp.Rational(1, den[i])*f5[i]
        D += c5
    if (depth >= 6):
        if debug == True:
            printd("BCH d'ordre 6")
        den = sp.S([-720, 240, 1440, 1440])
        f6 = [sp.S(0)]*len(den)
        f6[0] = Claudator(A, f5[2]) 
        f6[1] = Claudator(A, f5[4])
        f6[2] = Claudator(A, f5[5])
        f6[3] = Claudator(B, f5[0]) 
        c6 = sp.S(0)
        for i in range(0, len(den)):
            f6[i] = f6[i].expand()
            c6 += sp.Rational(1, den[i])*f6[i]
        D += c6
    if (depth >= 7):
        if debug == True:
            printd("BCH d'ordre 7")
        den = sp.S([30240, 5040, -10080, 10080, 1008, 5040, -7560, 3360, 10080, -10080, -1260, -1680, 3360, -3360, -2520, 7560, 10080, -30240])          
        f7 = [sp.S(0)]*len(den)
        f7[0]  = Claudator(A, Claudator(A, f5[0]))
        f7[1]  = Claudator(A, Claudator(A, f5[3]))
        f7[2]  = Claudator(A, f6[1])
        f7[3]  = Claudator(A, f6[3]) 
        f7[4]  = Claudator(A, Claudator(B, f5[1]))
        f7[5]  = Claudator(A, Claudator(B, f5[2])) 
        f7[6]  = Claudator(A, Claudator(B, f5[3]))
        f7[7]  = Claudator(A, Claudator(B, f5[4]))
        f7[8]  = Claudator(A, Claudator(B, f5[5]))
        f7[9]  = Claudator(B, Claudator(A, f5[0])) 
        f7[10] = Claudator(B, Claudator(A, f5[3]))
        f7[11] = Claudator(B, f6[1]) 
        f7[12] = Claudator(B, f6[3]) 
        f7[13] = Claudator(B, Claudator(B, f5[1])) 
        f7[14] = Claudator(B, Claudator(B, f5[2])) 
        f7[15] = Claudator(B, Claudator(B, f5[3]))
        f7[16] = Claudator(B, Claudator(B, f5[4]))
        f7[17] = Claudator(B, Claudator(B, f5[5]))
        c7 = sp.S(0)
        for i in range(0, len(den)):
            f7[i] = f7[i].expand()
            c7 += sp.Rational(1, den[i])*f7[i]
        D += c7
    if (depth >= 8):
        if debug == True:
            printd("BCH d'ordre 8")
        den = sp.S([sp.Rational(-24192, 5), 2520, 20160, 15120, -2016, -20160, 20160, -10080, -60480, -60480, 20160, -5040, 20160])   
        f8 = [sp.S(0)]*len(den)
        f8[0]  = Claudator(A, Claudator(A, f6[2]))
        f8[1]  = Claudator(A, f7[5])
        f8[2]  = Claudator(A, f7[8])
        f8[3]  = Claudator(A, Claudator(B, f6[0]))
        f8[4]  = Claudator(A, f7[11]) 
        f8[5]  = Claudator(A, Claudator(B, f6[2]))
        f8[6]  = Claudator(A, f7[13])
        f8[7]  = Claudator(A, f7[14])
        f8[8]  = Claudator(A, f7[17])
        f8[9]  = Claudator(B, f7[0])
        f8[10] = Claudator(B, Claudator(A, Claudator(A, f5[1])))
        f8[11] = Claudator(B, f7[1])
        f8[12] = Claudator(B, f7[9])
        c8 = sp.S(0)
        for i in range(0, len(den)):
            f8[i] = f8[i].expand()
            c8 += sp.Rational(1, den[i])*f8[i]
        D += c8
    if (depth >= 9):
        if debug == True:
            printd("BCH d'ordre 9")
        den = sp.S([-1209600, -302400, 113400, -40320, -120960, sp.Rational(604800, 13), 24192, 30240, -60480, -20160, sp.Rational(-604800, 11), -151200, -20160, -10080, 40320, 18900, -20160, -17280, -120960, -302400, 302400, 50400, 120960, 20160, -40320, 10080, 30240, 151200, 302400, 20160, -30240, 60480, sp.Rational(-604800, 13), -90720, 120960, 453600, 302400, 1209600])
        f9 = [sp.S(0)]*len(den)
        f9[0]  = Claudator(A, Claudator(A, f7[0]))
        f9[1]  = Claudator(A, Claudator(A, Claudator(A, Claudator(A, f5[1]))))
        f9[2]  = Claudator(A, Claudator(A, Claudator(A, f6[0])))
        f9[3]  = Claudator(A, Claudator(A, f7[2]))
        f9[4]  = Claudator(A, f8[0]) 
        f9[5]  = Claudator(A, f8[2])
        f9[6]  = Claudator(A, Claudator(A, Claudator(B, Claudator(A, f5[1]))))
        f9[7]  = Claudator(A, f8[3])
        f9[8]  = Claudator(A, Claudator(A, f7[10]))
        f9[9]  = Claudator(A, f8[5])
        f9[10] = Claudator(A, f8[8])
        f9[11] = Claudator(A, f8[9])
        f9[12] = Claudator(A, f8[10])
        f9[13] = Claudator(A, Claudator(B, f7[2]))
        f9[14] = Claudator(A, Claudator(B, Claudator(A, f6[2])))
        f9[15] = Claudator(A, Claudator(B, f7[8]))
        f9[16] = Claudator(A, Claudator(B, f7[11]))
        f9[17] = Claudator(A, Claudator(B, Claudator(B, f6[2])))
        f9[18] = Claudator(A, Claudator(B, f7[15]))
        f9[19] = Claudator(A, Claudator(B, f7[17]))
        f9[20] = Claudator(B, Claudator(A, f7[0]))
        f9[21] = Claudator(B, Claudator(A, Claudator(A, Claudator(A, f5[1]))))
        f9[22] = Claudator(B, Claudator(A, Claudator(A, f6[0])))
        f9[23] = Claudator(B, Claudator(A, f7[4]))
        f9[24] = Claudator(B, Claudator(A, f7[12]))
        f9[25] = Claudator(B, f8[6])
        f9[26] = Claudator(B, f8[7])
        f9[27] = Claudator(B, f8[8])
        f9[28] = Claudator(B, f8[9])
        f9[29] = Claudator(B, Claudator(B, f7[3]))
        f9[30] = Claudator(B, Claudator(B, f7[6]))
        f9[31] = Claudator(B, Claudator(B, f7[7]))
        f9[32] = Claudator(B, f8[12])
        f9[33] = Claudator(B, Claudator(B, Claudator(B, f6[0])))
        f9[34] = Claudator(B, Claudator(B, f7[12]))
        f9[35] = Claudator(B, Claudator(B, f7[15]))
        f9[36] = Claudator(B, Claudator(B, f7[16]))
        f9[37] = Claudator(B, Claudator(B, f7[17]))
        c9 = sp.S(0)
        for i in range(0, len(den)):
            c9 += sp.Rational(1, den[i])*f9[i]
        D += c9
    return D
