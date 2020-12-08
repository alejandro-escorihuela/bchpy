#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 07-12-2020
# alex
# test_Eel.py

import sympy as sp
from sympy.physics.quantum import Operator, Commutator

if __name__ == "__main__":
    A = Operator("A")
    B = Operator("B")
    e11 = A
    e12 = B
    e21 = Commutator(A,B)
    e31 = Commutator(A, Commutator(A, B))
    e32 = Commutator(B, Commutator(A, B))
    e41 = Commutator(A, Commutator(A, Commutator(A, B)))
    e42 = Commutator(B, Commutator(A, Commutator(A, B)))
    e43 = Commutator(B, Commutator(B, Commutator(B, A)))
    e51 = Commutator(A, Commutator(A, Commutator(A, Commutator(A, B))))
    e52 = Commutator(B, Commutator(A, Commutator(A, Commutator(A, B))))
    e53 = Commutator(A, Commutator(A, Commutator(B, Commutator(B, A))))
    e54 = Commutator(B, Commutator(B, Commutator(A, Commutator(A, B))))
    e55 = Commutator(A, Commutator(B, Commutator(B, Commutator(B, A))))
    e56 = Commutator(B, Commutator(B, Commutator(B, Commutator(B, B))))
    e61 = Commutator(A, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B)))))
    e62 = Commutator(B, Commutator(A, Commutator(A, Commutator(A, Commutator(A, B)))))
    e63 = Commutator(A, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B)))))
    e64 = Commutator(A, Commutator(B, Commutator(B, Commutator(A, Commutator(A, B)))))
    e65 = Commutator(B, Commutator(B, Commutator(A, Commutator(A, Commutator(A, B)))))
    e66 = Commutator(A, Commutator(A, Commutator(B, Commutator(B, Commutator(B, A)))))
    e67 = Commutator(B, Commutator(A, Commutator(B, Commutator(B, Commutator(B, A)))))
    e68 = Commutator(A, Commutator(B, Commutator(B, Commutator(B, Commutator(B, A)))))
    e69 = Commutator(B, Commutator(B, Commutator(B, Commutator(B, Commutator(B, A)))))
    
    r = sp.expand((Commutator(e21, e31) + e53 + e52).doit())
    print("[E21, E31] + E53 + E52 =", r)
    r = sp.expand((Commutator(e21, e32) + e54 + e55).doit())
    print("[E21, E32] + E54 + E55 =", r)    
    r = sp.expand((Commutator(e21, e41) - e63 + e62).doit())
    print("[E21, E41] - E63 + E62 =", r)
    r = sp.expand((Commutator(e21, e42) + (sp.Rational(1, 3))*(e65 + e66)).doit())
    print("[E21, E42] + (E65 + E66)/3 =", r)
    r = sp.expand((Commutator(e21, e43) - e68 + e67).doit())
    print("[E21, E43] - E68 + E67 =", r)
    r = sp.expand((Commutator(e31, e32) + e64 + (sp.Rational(1, 3))*(sp.S(2)*e66 - e65)).doit())
    print("[E31, E32] + E64 + (2*E66 - E65)/3 =", r)
    r = sp.expand((Commutator(e21, e41) - e63 + e62).doit())
    print("[E21, E41] - E63 + E62 =", r)
    r = sp.expand((Commutator(e11, e53) + sp.S(2)*e63 - e62).doit())
    print("[E11, E53] + 2*E63 - E62 =", r)
    r = sp.expand((Commutator(e12, e53) + e64 + (sp.Rational(1, 3))*(e65 + e66)).doit())
    print("[E12, E53] + E64 + (E65 + E66)/3 =", r)
