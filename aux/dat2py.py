#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 19-12-2022
# alex
# dat2py.py

import numpy as np

if __name__ == "__main__":
    with open("./bchLyndon20.dat") as fbch:
        linies = fbch.readlines()
        bch = np.zeros((111013 + 1, 6), dtype = float)
        i = 1
        for lin in linies:
            linia = lin.split(" ")
            j = 0
            if linia[0] != "\n":
                for k in range(len(linia)):
                    if linia[k] != "":
                        bch[i][j] = float(linia[k].replace("\n", ""))
                        j += 1
                i +=1
    bch[1][5], bch[2][5] = 1, 1
    for i in range(3, len(bch)):
        bch[i][5] = bch[int(bch[i][1])][5] + bch[int(bch[i][2])][5]
    print("elemBCH = []")
    for i in range(len(bch)):
        print("elemBCH.append([%d, %d, %d, %f, %f, %d])" % (bch[i][0], bch[i][1], bch[i][2], bch[i][3], bch[i][4], bch[i][5]))
