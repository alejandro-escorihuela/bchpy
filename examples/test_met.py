#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 24-04-2023
# alex
# geneqSS.py

import numpy as np
import sympy as sp
import time as tm
import sys
sys.path.insert(0, '../')
from bchpy import *
import relations as rl

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

def baseq(w):
    wq = []
    wq.append(w[0])
    wq.append(w[1])
    wq.append(w[2])
    wq.append(w[3])
    wq.append(w[4][:3])
    wq.append(w[5][2:5])
    wq.append(w[6][3:6])
    wq.append(w[7][5:11])
    wq.append(w[8][9:15])
    wq.append(w[9][17:27])
    return wq    



if __name__ == "__main__":
    tol = 1e-10

    ## Mètode d'escissió ABA complex sim-conj d'ordre 4
    cofs = [0.03, (0.089217853062385+0.025590308052868504j), 0.15, (0.18525002244352806-0.08754687281173089j), 0.20488095623525243, (0.22553212449408694+0.13887845131492216j),
            0.2302380875294951, (0.22553212449408694-0.13887845131492216j), 0.20488095623525243, (0.18525002244352806+0.08754687281173089j), 0.15,
            (0.089217853062385-0.025590308052868504j), 0.03]
    met = Metode(depth = 9, basis_type = "E", rkn = False)
    t0 = tm.time()
    met.setABA(*cofs)
    print("Temps ABA sim-conj d'ordre 4     = %f" % (tm.time() - t0))
    correcte = abs(met.w[1][1].real - 1.0) < tol and abs(met.w[1][1].imag) < tol
    correcte = correcte and abs(met.w[1][2].real - 1.0) < tol and abs(met.w[1][2].imag) < tol
    for i in range(2, 5):
        for j in range(1, len(met.w[i])):
            correcte = correcte and abs(met.w[i][j]) < tol
    if not correcte:
        printe("Resultat erroni en Mètode ABA complex sim-conj d'ordre 4")
        met.cprint()

    ## Mètode d'escissió BAB complex sim-conj d'ordre 4
    cofs = [0.03, (0.089217853062385+0.025590308052868504j), 0.15, (0.18525002244352806-0.08754687281173089j), 0.20488095623525243, (0.22553212449408694+0.13887845131492216j),
            0.2302380875294951, (0.22553212449408694-0.13887845131492216j), 0.20488095623525243, (0.18525002244352806+0.08754687281173089j), 0.15,
            (0.089217853062385-0.025590308052868504j), 0.03]
    met = Metode(depth = 9, basis_type = "E", rkn = False)
    t0 = tm.time()
    met.setBAB(*cofs)
    print("Temps BAB sim-conj d'ordre 4     = %f" % (tm.time() - t0))
    correcte = abs(met.w[1][1].real - 1.0) < tol and abs(met.w[1][1].imag) < tol
    correcte = correcte and abs(met.w[1][2].real - 1.0) < tol and abs(met.w[1][2].imag) < tol
    for i in range(2, 5):
        for j in range(1, len(met.w[i])):
            correcte = correcte and abs(met.w[i][j]) < tol
    if not correcte:
        printe("Resultat erroni en Mètode BAB complex sim-conj d'ordre 4")
        met.cprint()
        
    ## Mètode d'escissió RKN ABA d'ordre 8
    cofs = [0.052092434384033901862, 0.14585030481264474322, 0.22528749326770217132, 0.2551565441392939504,
            0.41627618961225709704,  0.018133468820831725316, -0.38456727021395037402, -0.17904011029926455989,
            0.09972717834705148443, -0.11847080143330224189, -0.108833834399100215506, 0.1864616892738210907,
            0.22201073664899168003, 0.45904158176713683037, 0.52387952203673426865, -0.0036608362703183590543,
            -0.5458724496837200138, -0.52694368162168635835, -0.5458724496837200138, -0.0036608362703183590543,
            0.52387952203673426865, 0.45904158176713683037, 0.22201073664899168003, 0.1864616892738210907,
            -0.108833834399100215506, -0.11847080143330224189, 0.09972717834705148443, -0.17904011029926455989,
            -0.38456727021395037402, 0.018133468820831725316, 0.41627618961225709704, 0.2551565441392939504,
            0.22528749326770217132, 0.14585030481264474322, 0.052092434384033901862]
    met = Metode(depth = 9, basis_type = "E", rkn = True)
    t0 = tm.time()
    met.setABA(*cofs)
    print("Temps RKN ABA sim d'ordre 8      = %f" % (tm.time() - t0))
    correcte = abs(met.w[1][1] + met.w[1][2] - 2.0) < tol
    wrkn = baserkn(met.w)
    for i in range(2, 9):
        for j in range(len(wrkn[i])):
            correcte = correcte and abs(wrkn[i][j]) < tol
    if not correcte:
        printe("Resultat erroni en Mètode RKN ABA simètric d'ordre 8")
    met.cprint()
    exit(-1)
    # ## Mètode d'escissió QA d'ordre 10
    # cofs = [0.3408726997217003, 0.6308398907950986, 0.5681165147290248, 0.6358118974865499, -0.194520172938841, -0.9046834271927621,
    #         -0.6985031964489967, -0.4814340856952135, 0.5611473127614152, 0.9689829571497421, 0.7120730850545042, 0.4317725809642063,
    #         0.00720904268022022, -0.40182836826909507, -0.5635668124920818, -0.800022967832949, -0.6820859994059998, 0.16950505672124352,
    #         0.74122378181587, 0.743992980440794, -0.3180809750107615, -0.49293651456761456, 0.052229439067892125]
    # cofs = cofs + cofs[-2::-1]
    # met = Metode(depth = 9, basis_type = "E", rkn = True)
    # t0 = tm.time()
    # met.setABA(*cofs)
    # print("Temps QA sim d'ordre 10          = %f" % (tm.time() - t0))
    # correcte = abs(met.w[1][1] + met.w[1][2] - 2.0) < tol
    # wq = baseq(met.w)
    # for i in range(2, 10):
    #     for j in range(len(wq[i])):
    #         correcte = correcte and abs(wq[i][j]) < tol
    # if not correcte:
    #     printe("Resultat erroni en Mètode QA simètric d'ordre 10")
    #     met.cprint()
        
    ## Mètode d'escissió ABA d'ordre 6
    cofs = [0.0502627644003922, 0.148816447901042, 0.413514300428344, -0.132385865767784, 0.0450798897943977, 0.067307604692185,
            -0.188054853819569, 0.432666402578175, 0.54196067845078, -0.016404589403617997, -0.7255255585086897, -0.016404589403617997,
            0.54196067845078, 0.432666402578175, -0.188054853819569, 0.067307604692185, 0.0450798897943977, -0.132385865767784,
            0.413514300428344, 0.148816447901042, 0.0502627644003922]
    met = Metode(depth = 9, basis_type = "E")
    t0 = tm.time()
    met.setABA(*cofs)
    print("Temps ABA simètric d'ordre 6     = %f" % (tm.time() - t0))
    correcte = abs(met.w[1][1] + met.w[1][2] - 2.0) < tol
    for i in range(2, 7):
        for j in range(1, len(met.w[i])):
            correcte = correcte and abs(met.w[i][j]) < tol
    if not correcte:
        printe("Resultat erroni en Mètode ABA simètric d'ordre 6")
        met.cprint()

    ## Mètode d'escissió BAB d'ordre 6
    cofs = [0.0502627644003922, 0.148816447901042, 0.413514300428344, -0.132385865767784, 0.0450798897943977, 0.067307604692185,
            -0.188054853819569, 0.432666402578175, 0.54196067845078, -0.016404589403617997, -0.7255255585086897,
            -0.016404589403617997, 0.54196067845078, 0.432666402578175, -0.188054853819569, 0.067307604692185,
            0.0450798897943977, -0.132385865767784, 0.413514300428344, 0.148816447901042, 0.0502627644003922]
    met = Metode(depth = 9, basis_type = "E")
    t0 = tm.time()
    met.setBAB(*cofs)
    print("Temps BAB simètric d'ordre 6     = %f" % (tm.time() - t0))
    correcte = abs(met.w[1][1] + met.w[1][2] - 2.0) < tol
    for i in range(2, 7):
        for j in range(1, len(met.w[i])):
            correcte = correcte and abs(met.w[i][j]) < tol
    if not correcte:
        printe("Resultat erroni en Mètode BAB simètric d'ordre 6")
        met.cprint()
    ## Mètode d'escissió SS->ABA d'ordre 10
    cofs = [0.039397861260843207443, 0.078795722521686414885, 0.19594591296839746625, 0.3130961034151085176,
            0.17050724332509329187, 0.027918383235078066129, -0.100837229179414508, -0.22959284159390708213,
            -0.049315390258371102905, 0.13096206107716487632, -0.069385672288672914365, -0.26973340565451070505,
            -0.09738003124930963428, 0.074973343155891436496, 0.093483383577850819035, 0.111993423999810201575,
            0.2390634367730184659, 0.36613344954622673022, -0.016486090294904581777, -0.39910563013603589377,
            -0.14800911580428240816, 0.10308739852747107746, 0.25725913624168065402, 0.4114308739558902306,
            0.20328225668637748451, -0.0048663605831352615624, -0.1984498571458875797, -0.39203335370863989784,
            -0.17004542537309512429, 0.051942502962449649262, 0.051303796861187071776, 0.05066509075992449429,
            0.050169730699827184972, 0.049674370639729875654, 0.049496053199662205097, 0.04931773575959453454,
            0.049496053199662205097, 0.049674370639729875654, 0.050169730699827184972, 0.05066509075992449429,
            0.051303796861187071776, 0.051942502962449649262, -0.17004542537309512429, -0.39203335370863989784,
            -0.1984498571458875797, -0.0048663605831352615624, 0.20328225668637748451, 0.4114308739558902306,
            0.25725913624168065402, 0.10308739852747107746, -0.14800911580428240816, -0.39910563013603589377,
            -0.016486090294904581777, 0.36613344954622673022, 0.2390634367730184659, 0.111993423999810201575,
            0.093483383577850819035, 0.074973343155891436496, -0.09738003124930963428, -0.26973340565451070505,
            -0.069385672288672914365, 0.13096206107716487632, -0.049315390258371102905, -0.22959284159390708213,
            -0.100837229179414508, 0.027918383235078066129, 0.17050724332509329187, 0.3130961034151085176,
            0.19594591296839746625, 0.078795722521686414885, 0.039397861260843207443]    
    met = Metode(depth = 9, basis_type = "E")
    t0 = tm.time()
    met.setABA(*cofs)
    print("Temps SS->ABA d'ordre 10         = %f" % (tm.time() - t0))
    correcte = abs(met.w[1][1] + met.w[1][2] - 2.0) < tol
    for i in range(2, 10):
        for j in range(1, len(met.w[i])):
            correcte = correcte and abs(met.w[i][j]) < tol
    if not correcte:
        printe("Resultat erroni en Mètode BAB simètric d'ordre 6")
        met.cprint()
        
    ## Mètode de composició XX d'ordre 6
    cofs = [0.0502627644003922, 0.0985536835006498, 0.31496061692769417, -0.44734648269547816, 0.49242637248987586,
            -0.42511876779769087, 0.23706391397812188, 0.19560248860005314, 0.34635818985072686, -0.36276277925434486,
            -0.36276277925434486, 0.34635818985072686, 0.19560248860005314, 0.23706391397812188, -0.42511876779769087,
            0.49242637248987586, -0.44734648269547816, 0.31496061692769417, 0.0985536835006498, 0.0502627644003922]
    met = Metode(depth = 9, basis_type = "C")
    t0 = tm.time()
    met.setXX(*cofs)
    print("Temps XX simètric d'ordre 6      = %f" % (tm.time() - t0))
    correcte = abs(met.w[1][1] - 1.0) < tol
    for i in range(2, 7):
        for j in range(1, len(met.w[i])):
            correcte = correcte and abs(met.w[i][j]) < tol
    if not correcte:
        printe("Resultat erroni en Mètode XX simètric d'ordre 6")
        met.cprint()
        
    ## Mètode SS(S2) simètric d'ordre 6
    cofs = [0.78451361047755726381949763, 0.23557321335935813368479318, -1.17767998417887100694641568, 1.31518632068391121888424973]
    met = Metode(depth = 13, basis_type = "M")
    t0 = tm.time()
    met.setSS2(cofs[0], cofs[1], cofs[2], cofs[3], cofs[2], cofs[1], cofs[0])
    print("Temps SS(S2) simètric d'ordre 6  = %f" % (tm.time() - t0))
    correcte = abs(met.w[1][1] - 1.0) < tol 
    for i in range(2, 7):
        for j in range(1, len(met.w[i])):
            correcte = correcte and abs(met.w[i][j]) < tol
    if not correcte:
        printe("Resultat erroni en Mètode SS(S2) simètric d'ordre 6")
        met.cprint()

    ## Mètode SS(S2) simètric d'ordre 8
    cofs = [0.74167036435061295344822780, -0.40910082580003159399730010, 0.19075471029623837995387626, -0.57386247111608226665638773,
            0.29906418130365592384446354, 0.33462491824529818378495798, 0.31529309239676659663205666, -0.79688793935291635291635401978884]
    cofs = cofs + cofs[-2::-1]
    met = Metode(depth = 13, basis_type = "M")
    t0 = tm.time()
    met.setSS2(*cofs)
    print("Temps SS(S2) simètric d'ordre 8  = %f" % (tm.time() - t0))
    correcte = abs(met.w[1][1] - 1.0) < tol
    for i in range(2, 9):
        for j in range(1, len(met.w[i])):
            correcte = correcte and abs(met.w[i][j]) < tol
    if not correcte:
        printe("Resultat erroni en Mètode SS(S2) simètric d'ordre 8")
        met.cprint()

    ## Mètode SS(S2) simètric d'ordre 10
    cofs = [0.078795722521686419263907679337684, 0.31309610341510852776481247192647, 0.027918383235078066109520273275299,
            -0.22959284159390709415121339679655, 0.13096206107716486317465685927961, -0.26973340565451071434460973222411,
            0.074973343155891435666137105641410, 0.11199342399981020488957508073640, 0.36613344954622675119314812353150,
            -0.39910563013603589787862981058340, 0.10308739852747107731580277001372, 0.41143087395589023782070411897608,
            -0.0048663605831352617621956593099771, -0.39203335370863990644808193642610, 0.051942502962449647037182904015976,
            0.050665090759924496335874344156866, 0.049674370639729879054568800279461, 0.049317735759594537917680008339338]
    cofs = cofs + cofs[-2::-1]
    met = Metode(depth = 13, basis_type = "M")
    t0 = tm.time()
    met.setSS2(*cofs)
    print("Temps SS(S2) simètric d'ordre 10 = %f" % (tm.time() - t0))
    correcte = abs(met.w[1][1] - 1.0) < tol
    for i in range(2, 10):
        for j in range(1, len(met.w[i])):
            correcte = correcte and abs(met.w[i][j]) < tol
    if not correcte:
        printe("Resultat erroni en Mètode SS(S2) simètric d'ordre 10")
        met.cprint()
        
    ## Mètode SC(S2) simètric-conjugat d'ordre 5
    cofs = [(0.1752684090720741140583563+0.05761474413053870201304364j), (0.1848736801929841604288898-0.1941219227572495885067758j),
            0.2797158214698834510255077,
            (0.1848736801929841604288898+0.1941219227572495885067758j), (0.1752684090720741140583563-0.05761474413053870201304364j)] 
    met = Metode(depth = 13, basis_type = "M")
    t0 = tm.time()
    met.setSS2(*cofs)
    print("Temps SC(S2) d'ordre 5           = %f" % (tm.time() - t0))
    correcte = abs(met.w[1][1].real - 1.0) < tol and abs(met.w[1][1].imag) < tol
    for i in range(2, 6):
        for j in range(1, len(met.w[i])):
            correcte = correcte and abs(met.w[i][j].real) < tol and abs(met.w[i][j].imag) < tol
    if not correcte:
        printe("Resultat erroni en Mètode SC(S2) simètric-conjugat d'ordre 5")
        met.cprint()        

    ## Mètode SS(S2) simètric d'ordre 6 amb coeficients complexos
    cofs =[(0.11690003755466129+0.04342825461606034j), (0.12955910128208825-0.1239896121880926j), (0.18653249281213383+0.003107430710072675j),
           (0.13401673670223335+0.15490785372391916j), (0.18653249281213383+0.003107430710072675j), (0.12955910128208825-0.1239896121880926j),
           (0.11690003755466129+0.04342825461606034j)]
    met = Metode(depth = 13, basis_type = "M")
    t0 = tm.time()
    met.setSS2(*cofs)
    print("Temps SS(S2) complex d'ordre 6   = %f" % (tm.time() - t0))
    correcte = abs(met.w[1][1].real - 1.0) < tol and abs(met.w[1][1].imag) < tol
    for i in range(2, 7):
        for j in range(1, len(met.w[i])):
            correcte = correcte and abs(met.w[i][j].real) < tol and abs(met.w[i][j].imag) < tol
    if not correcte:
        printe("Resultat erroni en Mètode SS(S2) simètric d'ordre 6 amb coeficients complexos")
        met.cprint()

    ## Mètode SS(S4) simètric d'ordre 6
    alpha = 1/(2 - 2**(1/5))
    beta = 1 -2*alpha
    cofs = [alpha, beta, alpha]
    met = Metode(depth = 15, basis_type = "M4")
    t0 = tm.time()
    met.setSS4(*cofs)
    print("Temps SS(S4) simètric d'ordre 6  = %f" % (tm.time() - t0))
    correcte = abs(met.w[1][1] - 1.0) < tol
    for i in range(2, 7):
        for j in range(1, len(met.w[i])):
            correcte = correcte and abs(met.w[i][j]) < tol
    if not correcte:
        printe("Resultat erroni en Mètode SS(S4) simètric d'ordre 6")
        met.cprint()
        
    ## Mètode SS(S4) simètric d'ordre 12
    cofs = [0.13699751300751958, 0.09957852698465418, 0.1800905109150247, 0.3431063349680623, -0.3875343134532472,
            -0.3067591226259879, 0.01556601702834869, 0.3217845977320384, 0.12953695882744667, -0.2951269901730427,
            0.06209067112774111, 0.4013385913228842, 0.06209067112774111, -0.2951269901730427, 0.12953695882744667,
            0.3217845977320384, 0.01556601702834869, -0.3067591226259879, -0.3875343134532472, 0.3431063349680623,
            0.1800905109150247, 0.09957852698465418, 0.13699751300751958]
    met = Metode(depth = 15, basis_type = "M4")
    t0 = tm.time()
    met.setSS4(*cofs)
    print("Temps SS(S4) d'ordre 4           = %f" % (tm.time() - t0))
    correcte = abs(met.w[1][1].real - 1.0) < tol and abs(met.w[1][1].imag) < tol
    for i in range(2, 13):
        for j in range(1, len(met.w[i])):
            correcte = correcte and abs(met.w[i][j].real) < tol and abs(met.w[i][j].imag) < tol
    if not correcte:
        printe("Resultat erroni en Mètode SS(S4) simètric d'ordre 12")
        met.cprint()
        
    ## Mètode SC(S4) d'ordre 12
    cofs = [(0.08037809771800615+0.0337165843080119j), (0.1093940318845407-0.01482087341961085j), (0.09262372956348475-0.06605710656363314j),
            (0.040995578476915186-0.11107926414727864j), (0.09969141287949691+0.0448722393657375j), (0.07626020627149126+0.09478617925994205j),
            0.0013138864121300880355,
            (0.07626020627149126-0.09478617925994205j), (0.09969141287949691-0.0448722393657375j), (0.040995578476915186+0.11107926414727864j),
            (0.09262372956348475+0.06605710656363314j), (0.1093940318845407+0.01482087341961085j), (0.08037809771800615-0.0337165843080119j)]
    met = Metode(depth = 15, basis_type = "M4")
    t0 = tm.time()
    met.setSS4(*cofs)
    print("Temps SC(S4) d'ordre 7(12)       = %f" % (tm.time() - t0))
    correcte = abs(met.w[1][1].real - 1.0) < tol and abs(met.w[1][1].imag) < tol
    for i in range(2, 8):
        for j in range(1, len(met.w[i])):
            correcte = correcte and abs(met.w[i][j].real) < tol and abs(met.w[i][j].imag) < tol
    if not correcte:
        printe("Resultat erroni en Mètode SC(S4) d'ordre 7(12)")
        met.cprint()

    ## Mètode SS(S4)->SS(S2) simètric d'ordre 12
    cofs = [0.18511202485634723, -0.23322653670517488, 0.18511202485634723, 0.13455122182641358, -0.16952391666817299, 0.13455122182641358, 0.24333959355207063,
            -0.30658867618911656, 0.24333959355207063, 0.4636077474157652, -0.5841091598634681, 0.4636077474157652, -0.523639151469176, 0.6597439894851047,
            -0.523639151469176, -0.4144951326914695, 0.522231142756951, -0.4144951326914695, 0.02103291415887125, -0.02649981128939381, 0.02103291415887125,
            0.43479766271737563, -0.5478107277027129, 0.43479766271737563, 0.17503127039222802, -0.2205255819570094, 0.17503127039222802, -0.39877771166322257,
            0.5024284331534025, -0.39877771166322257, 0.08389736138140566, -0.10570405163507021, 0.08389736138140566, 0.542291591006439, -0.6832445906899939,
            0.542291591006439, 0.08389736138140566, -0.10570405163507021, 0.08389736138140566, -0.39877771166322257, 0.5024284331534025, -0.39877771166322257,
            0.17503127039222802, -0.2205255819570094, 0.17503127039222802, 0.43479766271737563, -0.5478107277027129, 0.43479766271737563, 0.02103291415887125,
            -0.02649981128939381, 0.02103291415887125, -0.4144951326914695, 0.522231142756951, -0.4144951326914695, -0.523639151469176, 0.6597439894851047,
            -0.523639151469176, 0.4636077474157652, -0.5841091598634681, 0.4636077474157652, 0.24333959355207063, -0.30658867618911656, 0.24333959355207063,
            0.13455122182641358, -0.16952391666817299, 0.13455122182641358, 0.18511202485634723, -0.23322653670517488, 0.18511202485634723]
    met = Metode(depth = 13, basis_type = "M")
    t0 = tm.time()
    met.setSS2(*cofs)
    print("Temps SS(S4)->SS(S2) d'ordre 12  = %f" % (tm.time() - t0))
    correcte = abs(met.w[1][1].real - 1.0) < tol and abs(met.w[1][1].imag) < tol
    for i in range(2, 13):
        for j in range(1, len(met.w[i])):
            correcte = correcte and abs(met.w[i][j].real) < tol and abs(met.w[i][j].imag) < tol
    if not correcte:
        printe("Resultat erroni en Mètode SS(S4)->SS(S2) simètric d'ordre 12")
        met.cprint()
