#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 07-12-2020
# alex
# test_relE.py

import sympy as sp
from sympy.physics.quantum import Operator, Commutator
import sys
sys.path.insert(0, '../')
import relations as rl

if __name__ == "__main__":
    A = Operator("A")
    B = Operator("B")
    E = [[0 for j in range(rl.tamE[-1] + 1)] for i in range(len(rl.tamE) + 1)]
    E[1][1] = A
    E[1][2] = B
    E[2][1] = Commutator(A,B)
    E[3][1] = Commutator(A, Commutator(A, B))
    E[3][2] = Commutator(B, Commutator(A, B))
    E[4][1] = Commutator(A, Commutator(A, Commutator(A, B)))
    E[4][2] = Commutator(B, Commutator(A, Commutator(A, B)))
    E[4][3] = Commutator(B, Commutator(B, Commutator(B, A)))
    E[5][1] = Commutator(A, Commutator(A, Commutator(A, Commutator(A, B))))
    E[5][2] = Commutator(B, Commutator(A, Commutator(A, Commutator(A, B))))
    E[5][3] = Commutator(A, Commutator(A, Commutator(B, Commutator(B, A))))
    E[5][4] = Commutator(B, Commutator(B, Commutator(A, Commutator(A, B))))
    E[5][5] = Commutator(A, Commutator(B, Commutator(B, Commutator(B, A))))
    E[5][6] = Commutator(B, Commutator(B, Commutator(B, Commutator(B, A))))
    E[6][1] = Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B)))))
    E[6][2] = Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B)))))
    E[6][3] = Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B)))))
    E[6][4] = Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, B)))))
    E[6][5] = Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B)))))
    E[6][6] = Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(B, A)))))
    E[6][7] = Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(B, A)))))
    E[6][8] = Commutator(A, Commutator(B, Commutator(B, Commutator(B, Commutator(B, A)))))
    E[6][9] = Commutator(B, Commutator(B, Commutator(B, Commutator(B, Commutator(B, A)))))
    for i in range(1, 10):
        E[7][2*i - 1] = Commutator(A, E[6][i])  
        E[7][2*i] = Commutator(B, E[6][i])
    E[8][1]=Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B)))))))
    E[8][2]=Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B)))))))
    E[8][3]=Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B)))))))
    E[8][4]=Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B)))))))
    E[8][5]=Commutator(A, Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B)))))))
    E[8][6]=Commutator(B, Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B)))))))
    E[8][7]=Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B)))))))
    E[8][8]=Commutator(B, Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B)))))))
    E[8][9]=Commutator(B, Commutator(A, Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B)))))))
    E[8][10]=Commutator(A, Commutator(B, Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B)))))))
    E[8][11]=Commutator(B, Commutator(B, Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B)))))))
    E[8][12]=Commutator(A, Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, B)))))))
    E[8][13]=Commutator(B, Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, B)))))))
    E[8][14]=Commutator(A, Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, B)))))))
    E[8][15]=Commutator(B, Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, B)))))))
    E[8][16]=Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B)))))))
    E[8][17]=Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B)))))))
    E[8][18]=Commutator(A, Commutator(B, Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B)))))))
    E[8][19]=Commutator(B, Commutator(B, Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B)))))))
    E[8][20]=Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(B, A)))))))
    E[8][21]=Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(B, A)))))))
    E[8][22]=Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(B, A)))))))
    E[8][23]=Commutator(B, Commutator(A, Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(B, A)))))))
    E[8][24]=Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(B, A)))))))
    E[8][25]=Commutator(B, Commutator(B, Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(B, A)))))))
    E[8][26]=Commutator(B, Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(B, Commutator(B, A)))))))
    E[8][27]=Commutator(A, Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(B, Commutator(B, A)))))))
    E[8][28]=Commutator(B, Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(B, Commutator(B, A)))))))
    E[8][29]=Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(B, Commutator(B, Commutator(B, A)))))))
    E[8][30]=Commutator(B, Commutator(B, Commutator(B, Commutator(B, Commutator(B, Commutator(B, Commutator(B, A)))))))
    E[9][1]=Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][2]=Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][3]=Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][4]=Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][5]=Commutator(A, Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][6]=Commutator(B, Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][7]=Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][8]=Commutator(B, Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][9]=Commutator(A, Commutator(A, Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][10]=Commutator(B, Commutator(A, Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][11]=Commutator(A, Commutator(B, Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][12]=Commutator(B, Commutator(B, Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][13]=Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][14]=Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][15]=Commutator(A, Commutator(B, Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][16]=Commutator(B, Commutator(B, Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][17]=Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][18]=Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][19]=Commutator(A, Commutator(A, Commutator(B, Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][20]=Commutator(B, Commutator(A, Commutator(B, Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][21]=Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][22]=Commutator(B, Commutator(B, Commutator(B, Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][23]=Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, B))))))))
    E[9][24]=Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, B))))))))
    E[9][25]=Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, B))))))))
    E[9][26]=Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, B))))))))
    E[9][27]=Commutator(A, Commutator(A, Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, B))))))))
    E[9][28]=Commutator(B, Commutator(A, Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, B))))))))
    E[9][29]=Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, B))))))))
    E[9][30]=Commutator(B, Commutator(B, Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, B))))))))
    E[9][31]=Commutator(B, Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][32]=Commutator(A, Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][33]=Commutator(B, Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][34]=Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][35]=Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][36]=Commutator(A, Commutator(B, Commutator(B, Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][37]=Commutator(B, Commutator(B, Commutator(B, Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B))))))))
    E[9][38]=Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(B, A))))))))
    E[9][39]=Commutator(A, Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(B, A))))))))
    E[9][40]=Commutator(B, Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(B, A))))))))
    E[9][41]=Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(B, A))))))))
    E[9][42]=Commutator(B, Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(B, A))))))))
    E[9][43]=Commutator(A, Commutator(B, Commutator(A, Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(B, A))))))))
    E[9][44]=Commutator(B, Commutator(B, Commutator(A, Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(B, A))))))))
    E[9][45]=Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(B, A))))))))
    E[9][46]=Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(B, A))))))))
    E[9][47]=Commutator(A, Commutator(B, Commutator(B, Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(B, A))))))))
    E[9][48]=Commutator(B, Commutator(B, Commutator(B, Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(B, A))))))))
    E[9][49]=Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(B, Commutator(B, A))))))))
    E[9][50]=Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(B, Commutator(B, A))))))))
    E[9][51]=Commutator(B, Commutator(A, Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(B, Commutator(B, A))))))))
    E[9][52]=Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(B, Commutator(B, A))))))))
    E[9][53]=Commutator(B, Commutator(B, Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(B, Commutator(B, A))))))))
    E[9][54]=Commutator(B, Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(B, Commutator(B, Commutator(B, A))))))))
    E[9][55]=Commutator(A, Commutator(B, Commutator(B, Commutator(B, Commutator(B, Commutator(B, Commutator(B, Commutator(B, A))))))))
    E[9][56]=Commutator(B, Commutator(B, Commutator(B, Commutator(B, Commutator(B, Commutator(B, Commutator(B, Commutator(B, A))))))))
        
    for i in range(len(rl.tamE)):
        for j in range(rl.tamE[i - 1]):
            for k in range(len(rl.tamE)):
                for l in range(rl.tamE[k]):
                    if rl.relE[i + 1][j + 1][k + 1][l + 1][0][0] != 0:
                        res = rl.relE[i + 1][j + 1][k + 1][l + 1]
                        cad = "[E" + str(i + 1) + str(j + 1) + ", E" + str(k + 1) + str(l + 1) + "] - ("
                        primer = True
                        resq = sp.S(0)
                        for m in res:
                            resq += sp.Rational(m[0], m[1])*E[m[2]][m[3]]
                            if m[0] > 0 and primer == False:
                                cad += "+ "
                            elif m[0] < 0:
                                cad += "- "
                            primer = False
                            if abs(m[0]) > 1:
                                cad += str(abs(m[0])) + "*"
                            cad += "E" + str(m[2]) + str(m[3]) + ""
                            if m[1] != 1:
                                cad += "/" + str(m[1])
                            cad += " "
                        cad = cad[:-1] + ") ="
                        r = (Commutator(E[i + 1][j + 1], E[k + 1][l + 1]) - resq).doit().expand()
                        print(cad, r)
