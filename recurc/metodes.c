/* 29-01-2024 */
/* alex */
/* metodes.c */
#include <stdio.h>
#include <stdlib.h>
#include "metodes.h"

void rknitzarsim(double * res) {
  res[31] -= 0.6*res[37];
  res[33] += 1.6*res[37];
  res[34] += (2.0/15.0)*res[37];
  res[38] -= (4.0/3.0)*res[37];
  res[39] += 3*res[35] + res[58];
  res[40] -= res[37];
  res[41] -= 3*res[35] + 3*res[58];
  res[42] += res[35] + 0.5*res[58] - 1.5*res[59] - 3.5*res[61];
  res[43] -= (1.0/3.0)*res[55] + 0.5*res[60] + 0.5*res[62];
  res[44] += (4.0/3.0)*res[37];
  res[45] -= res[58];
  res[46] -= (4.0/3.0)*res[37];
  res[47] += 4*res[58];
  res[48] += 3*res[59] + 6*res[61];
  res[49] -= (5.0/9.0)*res[55] - (4.0/3.0)*res[60] - res[62];
  res[50] += res[37]/3;
  res[51] -= res[58];
  res[52] -= 3*res[59] + 9*res[61];
  res[53] += res[55];
  res[54] += 3*res[59] + 9*res[61];
  res[35] = res[37] = res[55] = res[58] = res[59] = res[60] = res[61] = res[62] = 0.0;
}

void rknitzar(double * res) {
  res[44] += 0.5*res[56];
  res[46] -= 0.5*res[56];
  res[47] -= 1.5*res[56];
  res[48] -= 1.5*res[57] + 3.5*res[58];
  res[49] -= res[56];
  res[50] += 4*res[56];
  res[51] += 3*res[57] + 6*res[58];
  res[52] -= res[56];
  res[53] -= 3*res[57] + 9*res[58];
  res[54] += 3*res[57] + 9*res[58];
  res[56] = res[57] = res[58] = 0.0;
  res[74] -= 0.6*res[80];
  res[76] += 1.6*res[80];
  res[77] += (2.0/15.0)*res[80];
  res[81] -= (4.0/3.0)*res[80];
  res[82] += 3*res[78] + res[101];
  res[83] -= res[80];
  res[84] -= 3*res[78] + 3*res[101];
  res[85] += res[78] + 0.5*res[101] - 1.5*res[102] - 3.5*res[104];
  res[86] -= (1.0/3.0)*res[98] + 0.5*res[103] + 0.5*res[105];
  res[87] += (4.0/3.0)*res[80];
  res[88] -= res[101];
  res[89] -= (4.0/3.0)*res[80];
  res[90] += 4*res[101];
  res[91] += 3*res[102] + 6*res[104];
  res[92] -= (5.0/9.0)*res[98] - (4.0/3.0)*res[103] - res[105];
  res[93] += res[80]/3;
  res[94] -= res[101];
  res[95] -= 3*res[102] + 9*res[104];
  res[96] += res[98];
  res[97] += 3*res[102] + 9*res[104];
  res[78] = res[80] = res[98] = res[101] = res[102] = res[103] = res[104] = res[105] = 0.0;
}

void metode_setABAsim(int tam, double * cofs, double * res, int order, int rkn) {
  int i, s = (tam - 1)/2;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0;
  if (s % 2 == 0)
    res[0] = cofs[s];
  else
    res[1] = cofs[s];
  for (i = s - 1; i >= 0; i--)
    if (i % 2 == 0)
      recAsim(res, cofs[i], order, rkn);
    else
      recBsim(res, cofs[i], order, rkn);
  if (rkn != 0)
    rknitzarsim(res);  
}

void metode_setBABsim(int tam, double * cofs, double * res, int order, int rkn) {
  int i, s = (tam - 1)/2;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0;
  if (s % 2 == 0)
    res[1] = cofs[s];
  else
    res[0] = cofs[s];
  for (i = s - 1; i >= 0; i--)
    if (i % 2 == 0)
      recBsim(res, cofs[i], order, rkn);
    else
      recAsim(res, cofs[i], order, rkn);
  if (rkn != 0)
    rknitzarsim(res);    
}

void metode_setABA(int tam, double * cofs, double * res, int order, int rkn) {
  int i;

  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0;
  res[0] = cofs[0];
  for (i = 1; i < tam; i++)
    if (i % 2 == 0)
      recA(res, cofs[i], order, rkn);
    else
      recB(res, cofs[i], order, rkn);
  if (rkn != 0)
    rknitzar(res);
}

void metode_setBAB(int tam, double * cofs, double * res, int order, int rkn) {
  int i;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0;
  res[1] = cofs[0];
  for (i = 1; i < tam; i++)
    if (i % 2 == 0)
      recB(res, cofs[i], order, rkn);
    else
      recA(res, cofs[i], order, rkn);
  if (rkn != 0)
    rknitzar(res);  
}

void metode_setXXsim(int tam, double * cofs, double * res, int order) {
  int i, s = tam/2;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0;
  for (i = s - 1; i >= 0; i--)
    if (i % 2 == 0)
      recXasim(res, cofs[i], order);
    else
      recXbsim(res, cofs[i], order);
}

void metode_setXX(int tam, double * cofs, double * res, int order) {
  int i;
  int vp[9] = {0, 1, 2, 4, 7, 13, 22, 40, 70};

  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0;
  res[vp[0]] = cofs[0];
  for (i = 1; i < order; i++)
    res[vp[i]] = res[vp[i - 1]]*cofs[0];
  for (i = 1; i < tam; i++)
    if (i % 2 == 0)
      recXa(res, cofs[i], order);
    else
      recXb(res, cofs[i], order); 
}

void metode_setS2sim(int tam, double * cofs, double * res, int order) {
  int i, m = tam/2;
  int vp[13] = {0, -1, 1, -1, 2, -1, 4, -1, 8, -1, 16, -1, 34};
  double auxp;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0;
  res[vp[0]] = cofs[m];
  auxp = cofs[m]*cofs[m];
  for (i = 2; i <= order; i += 2)
    res[vp[i]] = res[vp[i - 2]]*auxp;
  for (i = m - 1; i >= 0; i--)
    recS2sim(res, cofs[i], order);
}


void metode_setS2(int tam, double * cofs, double * res, int order) {
  int i;
  int vp[13] = {0, -1, 1, 2, 3, 5, 7, 11, 16, 24, 35, 53, 78};
  double auxp;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0;
  res[vp[0]] = cofs[0];
  auxp = cofs[0]*cofs[0];
  for (i = 2; i <= order; i += 2)
    res[vp[i]] = res[vp[i - 2]]*auxp;
  for (i = 1; i < tam; i++)
      recS2(res, cofs[i], order);  
}

void metode_setS4sim(int tam, double * cofs, double * res, int order) {
  int i, m = tam/2;
  int vp[15] = {0, -1, -1, -1, 1, -1, 2, -1, 4, -1, 7, -1, 12, -1, 21};
  double auxp;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0;
  res[vp[0]] = cofs[m];
  res[vp[4]] = pow(cofs[m], 5);
  auxp = cofs[m]*cofs[m];
  for (i = 6; i <= order; i+=2)
    res[vp[i]] = res[vp[i - 2]]*auxp;
  for (i = m - 1; i >= 0; i--)
      recS4sim(res, cofs[i], order);  
}


void metode_setS4(int tam, double * cofs, double * res, int order) {
  int i;
  int vp[15] = {0, -1, -1, -1, 1, 2, 3, 5, 7, 10, 13, 18, 24, 33, 44};
  double auxp;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0;
  res[vp[0]] = cofs[0];
  res[vp[4]] = pow(cofs[0], 5);
  auxp = cofs[0]*cofs[0];
  for (i = 6; i <= order; i+=2)
    res[vp[i]] = res[vp[i - 2]]*auxp;
  for (i = 1; i < tam; i++)
      recS4(res, cofs[i], order);  
}

void rknitzarsim_c(double complex * res) {
  res[31] -= 0.6*res[37];
  res[33] += 1.6*res[37];
  res[34] += (2.0/15.0)*res[37];
  res[38] -= (4.0/3.0)*res[37];
  res[39] += 3*res[35] + res[58];
  res[40] -= res[37];
  res[41] -= 3*res[35] + 3*res[58];
  res[42] += res[35] + 0.5*res[58] - 1.5*res[59] - 3.5*res[61];
  res[43] -= (1.0/3.0)*res[55] + 0.5*res[60] + 0.5*res[62];
  res[44] += (4.0/3.0)*res[37];
  res[45] -= res[58];
  res[46] -= (4.0/3.0)*res[37];
  res[47] += 4*res[58];
  res[48] += 3*res[59] + 6*res[61];
  res[49] -= (5.0/9.0)*res[55] - (4.0/3.0)*res[60] - res[62];
  res[50] += res[37]/3;
  res[51] -= res[58];
  res[52] -= 3*res[59] + 9*res[61];
  res[53] += res[55];
  res[54] += 3*res[59] + 9*res[61];
  res[35] = res[37] = res[55] = res[58] = res[59] = res[60] = res[61] = res[62] = 0.0 + 0.0*I;
}

void rknitzar_c(double complex * res) {
  res[44] += 0.5*res[56];
  res[46] -= 0.5*res[56];
  res[47] -= 1.5*res[56];
  res[48] -= 1.5*res[57] + 3.5*res[58];
  res[49] -= res[56];
  res[50] += 4*res[56];
  res[51] += 3*res[57] + 6*res[58];
  res[52] -= res[56];
  res[53] -= 3*res[57] + 9*res[58];
  res[54] += 3*res[57] + 9*res[58];
  res[56] = res[57] = res[58] = 0.0 + 0.0*I;
  res[74] -= 0.6*res[80];
  res[76] += 1.6*res[80];
  res[77] += (2.0/15.0)*res[80];
  res[81] -= (4.0/3.0)*res[80];
  res[82] += 3*res[78] + res[101];
  res[83] -= res[80];
  res[84] -= 3*res[78] + 3*res[101];
  res[85] += res[78] + 0.5*res[101] - 1.5*res[102] - 3.5*res[104];
  res[86] -= (1.0/3.0)*res[98] + 0.5*res[103] + 0.5*res[105];
  res[87] += (4.0/3.0)*res[80];
  res[88] -= res[101];
  res[89] -= (4.0/3.0)*res[80];
  res[90] += 4*res[101];
  res[91] += 3*res[102] + 6*res[104];
  res[92] -= (5.0/9.0)*res[98] - (4.0/3.0)*res[103] - res[105];
  res[93] += res[80]/3;
  res[94] -= res[101];
  res[95] -= 3*res[102] + 9*res[104];
  res[96] += res[98];
  res[97] += 3*res[102] + 9*res[104];
  res[78] = res[80] = res[98] = res[101] = res[102] = res[103] = res[104] = res[105] = 0.0 + 0.0*I;  
}

void metode_setABAsim_c(int tam, double complex * cofs, double complex * res, int order, int rkn) {
  int i, s = (tam - 1)/2;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0 + 0.0*I;
  if (s % 2 == 0)
    res[0] = cofs[s];
  else
    res[1] = cofs[s];
  for (i = s - 1; i >= 0; i--)
    if (i % 2 == 0)
      recAsim_c(res, cofs[i], order, rkn);
    else
      recBsim_c(res, cofs[i], order, rkn);
  if (rkn != 0)
    rknitzarsim_c(res);    
}

void metode_setBABsim_c(int tam, double complex * cofs, double complex * res, int order, int rkn) {
  int i, s = (tam - 1)/2;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0 + 0.0*I;
  if (s % 2 == 0)
    res[1] = cofs[s];
  else
    res[0] = cofs[s];
  for (i = s - 1; i >= 0; i--)
    if (i % 2 == 0)
      recBsim_c(res, cofs[i], order, rkn);
    else
      recAsim_c(res, cofs[i], order, rkn);
  if (rkn != 0)
    rknitzarsim_c(res);    
}

void metode_setABA_c(int tam, double complex * cofs, double complex * res, int order, int rkn) {
  int i;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0 + 0.0*I;
  res[0] = cofs[0];
  for (i = 1; i <= tam; i++)
    if (i % 2 == 0)
      recA_c(res, cofs[i], order, rkn);
    else
      recB_c(res, cofs[i], order, rkn);
  if (rkn != 0)
    rknitzar_c(res);  
}

void metode_setBAB_c(int tam, double complex * cofs, double complex * res, int order, int rkn) {
  int i;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0 + 0.0*I;
  res[1] = cofs[0];
  for (i = 1; i < tam; i++)
    if (i % 2 == 0)
      recB_c(res, cofs[i], order, rkn);
    else
      recA_c(res, cofs[i], order, rkn);
  if (rkn != 0)
    rknitzar_c(res);  
}

void metode_setXXsim_c(int tam, double complex * cofs, double complex * res, int order) {
  int i, s = tam/2;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0 + 0.0*I;
  for (i = s - 1; i >= 0; i--)
    if (i % 2 == 0)
      recXasim_c(res, cofs[i], order);
    else
      recXbsim_c(res, cofs[i], order);
}

void metode_setXX_c(int tam, double complex * cofs, double complex * res, int order) {
  int i;
  int vp[9] = {0, 1, 2, 4, 7, 13, 22, 40, 70};
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0 + 0.0*I;
  res[vp[0]] = cofs[0];
  for (i = 1; i < order; i++)
    res[vp[i]] = res[vp[i - 1]]*cofs[0];
  for (i = 1; i < tam; i++)
    if (i % 2 == 0)
      recXa_c(res, cofs[i], order);
    else
      recXb_c(res, cofs[i], order); 
}

void metode_setS2sim_c(int tam, double complex * cofs, double complex * res, int order) {
  int i, m = tam/2;
  int vp[13] = {0, -1, 1, -1, 2, -1, 4, -1, 8, -1, 16, -1, 34};
  double complex auxp;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0 + 0.0*I;
  res[vp[0]] = cofs[m];
  auxp = cofs[m]*cofs[m];
  for (i = 2; i <= order; i += 2)
    res[vp[i]] = res[vp[i - 2]]*auxp;
  for (i = m - 1; i >= 0; i--)
    recS2sim_c(res, cofs[i], order);
}

void metode_setS2_c(int tam, double complex * cofs, double complex * res, int order) {
  int i;
  int vp[13] = {0, -1, 1, 2, 3, 5, 7, 11, 16, 24, 35, 53, 78};
  double complex auxp;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0 + 0.0*I;
  res[vp[0]] = cofs[0];
  auxp = cofs[0]*cofs[0];
  for (i = 2; i <= order; i+=2)
    res[vp[i]] = res[vp[i - 2]]*auxp;
  for (i = 1; i < tam; i++)
      recS2_c(res, cofs[i], order);  
}

void metode_setS4sim_c(int tam, double complex * cofs, double complex * res, int order) {
  int i, m = tam/2;
  int vp[15] = {0, -1, -1, -1, 1, -1, 2, -1, 4, -1, 7, -1, 12, -1, 21};
  double complex auxp;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0;
  res[vp[0]] = cofs[m];
  res[vp[4]] = cpow(cofs[m], 5);
  auxp = cofs[m]*cofs[m];
  for (i = 6; i <= order; i+=2)
    res[vp[i]] = res[vp[i - 2]]*auxp;
  for (i = m - 1; i >= 0; i--)
      recS4sim_c(res, cofs[i], order);    
}


void metode_setS4_c(int tam, double complex * cofs, double complex * res, int order) {
  int i;
  int vp[15] = {0, -1, -1, -1, 1, 2, 3, 5, 7, 10, 13, 18, 24, 33, 44};
  double complex auxp;
  
  for (i = 0; i < TAMBCH; i++)
    res[i] = 0.0 + 0.0*I;
  res[vp[0]] = cofs[0];
  res[vp[4]] = cpow(cofs[0], 5);
  auxp = cofs[0]*cofs[0];
  for (i = 6; i <= order; i+=2)
    res[vp[i]] = res[vp[i - 2]]*auxp;
  for (i = 1; i < tam; i++)
      recS4_c(res, cofs[i], order);  
}
