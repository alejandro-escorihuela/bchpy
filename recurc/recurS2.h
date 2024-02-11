/* 05-02-2024 */
/* alex */
/* recurS2.h */
#ifndef _RECURS2_H
#define _RECURS2_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "esq.h"

void recS2sim(double * res, double x, int order);
void recS2(double * res, double x, int order);

/* Composició SS(S2)
Ordre 1
res[0] <- M_1,1

Ordre 3
res[1] <- M_3,1

Ordre 4
res[2] <- M_4,1

Ordre 5
res[3] <- M_5,1
res[4] <- M_5,2

Ordre 6
res[5] <- M_6,1
res[6] <- M_6,2

Ordre 7
res[7] <- M_7,1
res[8] <- M_7,2
res[9] <- M_7,3
res[10] <- M_7,4

Ordre 8
res[11] <- M_8,1
res[12] <- M_8,2
res[13] <- M_8,3
res[14] <- M_8,4
res[15] <- M_8,5

Ordre 9
res[16] <- M_9,1
res[17] <- M_9,2
res[18] <- M_9,3
res[19] <- M_9,4
res[20] <- M_9,5
res[21] <- M_9,6
res[22] <- M_9,7
res[23] <- M_9,8

Ordre 10
res[24] <- M_10,1
res[25] <- M_10,2
res[26] <- M_10,3
res[27] <- M_10,4
res[28] <- M_10,5
res[29] <- M_10,6
res[30] <- M_10,7
res[31] <- M_10,8
res[32] <- M_10,9
res[33] <- M_10,10
res[34] <- M_10,11

Ordre 11
res[35] <- M_11,1
res[36] <- M_11,2
res[37] <- M_11,3
res[38] <- M_11,4
res[39] <- M_11,5
res[40] <- M_11,6
res[41] <- M_11,7
res[42] <- M_11,8
res[43] <- M_11,9
res[44] <- M_11,10
res[45] <- M_11,11
res[46] <- M_11,12
res[47] <- M_11,13
res[48] <- M_11,14
res[49] <- M_11,15
res[50] <- M_11,16
res[51] <- M_11,17
res[52] <- M_11,18

Ordre 12
res[53] <- M_11,1
res[54] <- M_11,2
res[55] <- M_11,3
res[56] <- M_11,4
res[57] <- M_11,5
res[58] <- M_11,6
res[59] <- M_11,7
res[60] <- M_11,8
res[61] <- M_11,9
res[62] <- M_11,10
res[63] <- M_11,11
res[64] <- M_11,12
res[65] <- M_11,13
res[66] <- M_11,14
res[67] <- M_11,15
res[68] <- M_11,16
res[69] <- M_11,17
res[70] <- M_11,18
res[71] <- M_11,19
res[72] <- M_11,20
res[73] <- M_11,21
res[74] <- M_11,22
res[75] <- M_11,23
res[76] <- M_11,24
res[77] <- M_11,25

Ordre 13
res[78] <- M_11,1
res[79] <- M_11,2
res[80] <- M_11,3
res[81] <- M_11,4
res[82] <- M_11,5
res[83] <- M_11,6
res[84] <- M_11,7
res[85] <- M_11,8
res[86] <- M_11,9
res[87] <- M_11,10
res[88] <- M_11,11
res[89] <- M_11,12
res[90] <- M_11,13
res[91] <- M_11,14
res[92] <- M_11,15
res[93] <- M_11,16
res[94] <- M_11,17
res[95] <- M_11,18
res[96] <- M_11,19
res[97] <- M_11,20
res[98] <- M_11,21
res[99] <- M_11,22
res[100] <- M_11,23
res[101] <- M_11,24
res[102] <- M_11,25
res[103] <- M_11,26
res[104] <- M_11,27
res[105] <- M_11,28
res[106] <- M_11,29
res[107] <- M_11,30
res[108] <- M_11,31
res[109] <- M_11,32
res[110] <- M_11,33
res[111] <- M_11,34
res[112] <- M_11,35
res[113] <- M_11,36
res[114] <- M_11,37
res[115] <- M_11,38
res[116] <- M_11,39
res[117] <- M_11,40

*/
#endif
