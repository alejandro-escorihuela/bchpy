#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 01-01-2021
# alex
# qolib.py

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
        tot_num = True
        for elem in sol:
            tot_num = tot_num and sol[elem].is_number
        if tot_num == False:
            cpe_arg = cpe.args
            pos_afg = []
            for i in range(len(cpe_arg)):
                cnc_i = cpe_arg[i].args_cnc()
                equ = sp.S(0)
                for j in range(i, len(cpe_arg)):
                    cnc_j = cpe_arg[j].args_cnc()
                    if cnc_i[1] == cnc_j[1] and j not in pos_afg:
                        pos_afg.append(j)
                        equ += sp.prod(cnc_j[0])
                if (len(equ.args) > 0) and equ not in eqs and -equ not in eqs:
                    eqs.append(equ)
            sol = sp.solve(eqs)
    return sol

def escl(can, En):
    if len(En) == 1:
        return False
    sim = sp.symbols("s1:" + str(len(En)))
    cesq = sp.S(0)
    numA = can.count(A)
    for i in range(1, len(En)):
        if numA == En[i].count(A):
            cesq += sim[i - 1]*En[i]
    if cesq == 0:
        return False
    cesq = cesq.doit().expand()
    cane = can.doit().expand()
    if resoldre(cesq, cane):
        return True
    return False

def pcorxets(E):
    for i in range(1, len(E)):
        for j in range(1, len(E[i])):
            print("E_{" + str(i) + ", " + str(j) + "} =", E[i][j])

