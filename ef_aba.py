#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 28-06-2022
# alex
# ef.py

from bchpy import *

def baserkn(w):
    wrkn = []
    wrkn.append(w[0])
    wrkn.append(w[1])
    wrkn.append(w[2])
    wrkn.append(w[3])
    wrkn.append(w[4][:3])
    wrkn.append(w[5][:5])
    wrkn.append(w[6][:6])
    wrkn.append(w[7][:11])
    wrkn.append(w[8][:15])
    wrkn.append(w[9][1:10] + w[9][11:27])
    return wrkn

def ef(s, m, vec, sc):
    err = 0.0
    if sc:
        if m % 2 == 0:
            for i in range(len(vec)):
                err += np.imag(vec[i])**2
        else:
            for i in range(len(vec)):
                err += np.real(vec[i])**2            
    else:
        for i in range(len(vec)):
            err += vec[i]**2
    err = np.sqrt(float(err))
    err = err**(1/(m - 1))
    if err < 0.01:
        return 0.0
    else:
        return s*err

def eferr(cofs, ABA = True, sc = False, rkn = False):
    a, b = cofs
    met = Metode(depth = 9, basis_type = "E", numeric = True)
    cc = []
    if ABA:
        s = len(b)
        for i in range(s):
            cc.append(a[i])
            cc.append(b[i])
        cc.append(a[s])
        met.setABA(*cc)
    else:
        s = len(a)
        for i in range(s):
            cc.append(b[i])
            cc.append(a[i])
        cc.append(b[s])
        met.setBAB(*cc)        
    ret = []
    if rkn:
        base = baserkn(met.w)
    else:
        base = met.w
    for i in range(2, 10):
        ret.append(ef(s, i, base[i], sc))
    return ret
