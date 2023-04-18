#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 01-01-2021
# alex
# qolib.py

import sympy as sp
from sympy.physics.quantum import Operator, Commutator
from sympy.parsing.sympy_parser import parse_expr

def expOp(A, n):
    fac = sp.S(1)
    ex = sp.S(1)
    for i in range(1, n):
        fac *= sp.S(i)
        ex += A**i/fac
    return ex

def logOp(O, n):
    lg = sp.S(1)
    sig = -1
    for i in range(1, n):
        sig *= -1
        f = sp.Rational(sig, i)
        op = O - sp.S(1)
        op = delterms(op, Operator("A"), Operator("B"), n - i + 1)
        lg += f*op**i
    return lg

def bchOp(A, B, n):
    expo = sp.expand(expOp(A, n)*expOp(B, n))
    expo = delterms(expo, A, B, n)
    loga = sp.expand(logOp(expo, n))
    loga = delterms(loga, A, B, n)
    return loga

def countDim(expre, opA, opB):
    dim = 0
    acnc = expre.args_cnc()[1]
    for i in range(0, len(acnc)):
        if len(acnc[i].args) > 1:
            dim += acnc[i].args[1]
        else:
            dim += 1
    return dim

def delterms(expre, opA, opB, n):
    retex = sp.S(0)
    for arg in expre.args:
        if countDim(arg, opA, opB) < n:
            retex += arg
    return retex

def llegir_base(txt, opA, opB):
    f = open(txt, "r")
    linies = f.readlines()
    E = []
    E.append([0])
    i = 1
    for lin in linies:
        E.append([0])
        ls = lin.replace("\n", "").split(" ")[:-1]
        for j in range(0, len(ls)):
            item = ls[j].replace("[", "").replace("]", "")
            item_rev = item[::-1]
            com = opA
            if item_rev[0] == "B":
                com = opB
            for k in range(1, len(item_rev)):
                if item_rev[k] == "A":
                    com = Commutator(opA, com)
                else:
                    com = Commutator(opB, com)
            E[i].append(com)
        i += 1
    return E

def llegir_baseCAdj(txt, Y):
    f = open(txt, "r")
    linies = f.readlines()
    Z = []
    Z.append([0])
    i = 1
    for lin in linies:
        Z.append([0])
        ls = lin.replace("\n", "").split(" ")[:-1]
        for j in range(0, len(ls)):
            item = ls[j].replace("[", "Commutator(").replace("]", ")").replace("_{", "[").replace("}", "]")
            com = eval(item)
            Z[i].append(com)
        i += 1
    return Z

def llegir_baseS2(txt, S):
    f = open(txt, "r")
    linies = f.readlines()
    Z = []
    Z.append([0])
    i = 1
    for lin in linies:
        Z.append([0])
        ls = lin.replace("\n", "").split(" ")[:-1]
        if len(ls) == 0:
            Z[i].append(0)
        for j in range(0, len(ls)):
            item = ls[j].replace("[", "Commutator(").replace("]", ")").replace("_{", "[").replace("}", "]")
            com = eval(item)
            Z[i].append(com)
        i += 1
    return Z

def resoldre(exprEsq, exprDre):
    eqs = []
    cpe = exprEsq
    eDre_arg = exprDre.args
    eEsq_arg = exprEsq.args
    for i in range(len(eDre_arg)):
        cnc_i = eDre_arg[i].args_cnc()
        equ = sp.S(0)
        for j in range(len(eEsq_arg)):
            cnc_j = eEsq_arg[j].args_cnc()
            if cnc_i[1] == cnc_j[1]:
                equ += sp.prod(cnc_j[0])
                cpe -= eEsq_arg[j]
        equ -= sp.prod(cnc_i[0])
        if equ not in eqs and -equ not in eqs:
            eqs.append(equ)
    sol = sp.solve(eqs)
    if sol:
        cpe_arg = cpe.args
        pos_afg = []
        #print(cpe_arg)
        for i in range(len(cpe_arg)):
            cnc_i = cpe_arg[i].args_cnc()
            equ = sp.S(0)
            for j in range(i, len(cpe_arg)):
                cnc_j = cpe_arg[j].args_cnc()
                if cnc_i[1] == cnc_j[1] and j not in pos_afg:
                    pos_afg.append(j)
                    #print(cnc_j[0], cnc_j[1], sp.prod(cnc_j[0]))
                    equ += sp.prod(cnc_j[0])
            #print(equ, len(equ.args))
            #if (len(equ.args) > 0) and equ not in eqs and -equ not in eqs:
            if equ not in eqs and -equ not in eqs:
                eqs.append(equ)
        sol = sp.solve(eqs)
        #print(eqs)
    return sol

def esclCA(can, En, A):
    if len(En) == 1:
        return []
    sim = sp.symbols("s1:" + str(len(En)))
    cesq = sp.S(0)
    numA = can.count(A)
    for i in range(1, len(En)):
        if numA == En[i].count(A):
            cesq += sim[i - 1]*En[i]
    if cesq == 0:
        return []
    cesq = cesq.doit().expand()
    cane = can.doit().expand()
    return resoldre(cesq, cane)

def escl(can, En):
    if len(En) == 1:
        return []
    sim = sp.symbols("s1:" + str(len(En)))
    cesq = sp.S(0)
    for i in range(1, len(En)):
        cesq += sim[i - 1]*En[i]
    if cesq == 0:
        return []
    cesq = cesq.doit().expand()
    cane = can.doit().expand()
    return resoldre(cesq, cane)

def pcorxets(E, label = "N"):
    for i in range(1, len(E)):
        for j in range(1, len(E[i])):
            print(label + "_{" + str(i) + ", " + str(j) + "} =", E[i][j])

