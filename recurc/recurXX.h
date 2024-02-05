/* 03-02-2024 */
/* alex */
/* recurXX.h */
#ifndef _RECURXX_H
#define _RECURXX_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "esq.h"

void recXa(double * res, double x, int order);
void recXb(double * res, double x, int order);

/* Composició mètode i adjunt
Ordre 1
res[0] <- C_1,1

Ordre 2
res[1] <- C_2,1

Ordre 3
res[2] <- C_3,1
res[3] <- C_3,2

Ordre 4
res[4] <- C_4,1
res[5] <- C_4,2
res[6] <- C_4,3

Ordre 5
res[7] <- C_5,1
res[8] <- C_5,2
res[9] <- C_5,3
res[10] <- C_5,4
res[11] <- C_5,5
res[12] <- C_5,6

Ordre 6
res[13] <- C_6,1
res[14] <- C_6,2
res[15] <- C_6,3
res[16] <- C_6,4
res[17] <- C_6,5
res[18] <- C_6,6
res[19] <- C_6,7
res[20] <- C_6,8
res[21] <- C_6,9

Ordre 7
res[22] <- C_7,1
res[23] <- C_7,2
res[24] <- C_7,3
res[25] <- C_7,4
res[26] <- C_7,5
res[27] <- C_7,6
res[28] <- C_7,7
res[29] <- C_7,8
res[30] <- C_7,9
res[31] <- C_7,10
res[32] <- C_7,11
res[33] <- C_7,12
res[34] <- C_7,13
res[35] <- C_7,14
res[36] <- C_7,15
res[37] <- C_7,16
res[38] <- C_7,17
res[39] <- C_7,18

Ordre 8
res[40] <- C_8,1
res[41] <- C_8,2
res[42] <- C_8,3
res[43] <- C_8,4
res[44] <- C_8,5
res[45] <- C_8,6
res[46] <- C_8,7
res[47] <- C_8,8
res[48] <- C_8,9
res[49] <- C_8,10
res[50] <- C_8,11
res[51] <- C_8,12
res[52] <- C_8,13
res[53] <- C_8,14
res[54] <- C_8,15
res[55] <- C_8,16
res[56] <- C_8,17
res[57] <- C_8,18
res[58] <- C_8,19
res[59] <- C_8,20
res[60] <- C_8,21
res[61] <- C_8,22
res[62] <- C_8,23
res[63] <- C_8,24
res[64] <- C_8,25
res[65] <- C_8,26
res[66] <- C_8,27
res[67] <- C_8,28
res[68] <- C_8,29
res[69] <- C_8,30

Ordre 9
res[70] <- C_9,1
res[71] <- C_9,2
res[72] <- C_9,3
res[73] <- C_9,4
res[74] <- C_9,5
res[75] <- C_9,6
res[76] <- C_9,7
res[77] <- C_9,8
res[78] <- C_9,9
res[79] <- C_9,10
res[80] <- C_9,11
res[81] <- C_9,12
res[82] <- C_9,13
res[83] <- C_9,14
res[84] <- C_9,15
res[85] <- C_9,16
res[86] <- C_9,17
res[87] <- C_9,18
res[88] <- C_9,19
res[89] <- C_9,20
res[90] <- C_9,21
res[91] <- C_9,22
res[92] <- C_9,23
res[93] <- C_9,24
res[94] <- C_9,25
res[95] <- C_9,26
res[96] <- C_9,27
res[97] <- C_9,28
res[98] <- C_9,29
res[99] <- C_9,30
res[100] <- C_9,31
res[101] <- C_9,32
res[102] <- C_9,33
res[103] <- C_9,34
res[104] <- C_9,35
res[105] <- C_9,36
res[106] <- C_9,37
res[107] <- C_9,38
res[108] <- C_9,39
res[109] <- C_9,40
res[110] <- C_9,41
res[111] <- C_9,42
res[112] <- C_9,43
res[113] <- C_9,44
res[114] <- C_9,45
res[115] <- C_9,46
res[116] <- C_9,47
res[117] <- C_9,48
res[118] <- C_9,49
res[119] <- C_9,50
res[120] <- C_9,51
res[121] <- C_9,52
res[122] <- C_9,53
res[123] <- C_9,54
res[124] <- C_9,55
res[125] <- C_9,56
*/
  
#endif

