/* 28-02-2022 */
/* alex */
/* recur.h */
#ifndef _RECUR_H
#define _RECUR_H
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#define TAMBCH 128

void copyesq(double * dest, double * orig);

void recAsim_c(double * res, double x, int order, int rkn);
void recBsim_c(double * res, double x, int order, int rkn);

/*
Escissió no RKN

Ordre 1
res[0] <- E_1,1
res[1] <- E_1,2

Ordre 3
res[2] <- E_3,1
res[3] <- E_3,2

Ordre 5
res[4] <- E_5,1
res[5] <- E_5,2
res[6] <- E_5,3
res[7] <- E_5,4
res[8] <- E_5,5
res[9] <- E_5,6

Ordre 7
res[10] <- E_7,1
res[11] <- E_7,2
res[12] <- E_7,3
res[13] <- E_7,4
res[14] <- E_7,5
res[15] <- E_7,6
res[16] <- E_7,7
res[17] <- E_7,8
res[18] <- E_7,9
res[19] <- E_7,10
res[20] <- E_7,11
res[21] <- E_7,12
res[22] <- E_7,13
res[23] <- E_7,14
res[24] <- E_7,15
res[25] <- E_7,16
res[26] <- E_7,17
res[27] <- E_7,18

Ordre 9
res[28] <- E_9,1
res[29] <- E_9,2
res[30] <- E_9,3
res[31] <- E_9,4
res[32] <- E_9,5
res[33] <- E_9,6
res[34] <- E_9,7
res[35] <- E_9,8
res[36] <- E_9,9
res[37] <- E_9,10
res[38] <- E_9,11
res[39] <- E_9,12
res[40] <- E_9,13
res[41] <- E_9,14
res[42] <- E_9,15
res[43] <- E_9,16
res[44] <- E_9,17
res[45] <- E_9,18
res[46] <- E_9,19
res[47] <- E_9,20
res[48] <- E_9,21
res[49] <- E_9,22
res[50] <- E_9,23
res[51] <- E_9,24
res[52] <- E_9,25
res[53] <- E_9,26
res[54] <- E_9,27
res[55] <- E_9,28
res[56] <- E_9,29
res[57] <- E_9,30
res[58] <- E_9,31
res[59] <- E_9,32
res[60] <- E_9,33
res[61] <- E_9,34
res[62] <- E_9,35
res[63] <- E_9,36
res[64] <- E_9,37
res[65] <- E_9,38
res[66] <- E_9,39
res[67] <- E_9,40
res[68] <- E_9,41
res[69] <- E_9,42
res[70] <- E_9,43
res[71] <- E_9,44
res[72] <- E_9,45
res[73] <- E_9,46
res[74] <- E_9,47
res[75] <- E_9,48
res[76] <- E_9,49
res[77] <- E_9,50
res[78] <- E_9,51
res[79] <- E_9,52
res[80] <- E_9,53
res[81] <- E_9,54
res[82] <- E_9,55
res[83] <- E_9,56
*/

/*
Escissió RKN

Ordre 1
res[0] <- E_1,1
res[1] <- E_1,2

Ordre 3
res[2] <- E_3,1
res[3] <- E_3,2

Ordre 5
res[4] <- E_5,1
res[5] <- E_5,2
res[6] <- E_5,3
res[7] <- E_5,4

Ordre 7
res[10] <- E_7,1
res[11] <- E_7,2
res[12] <- E_7,3
res[13] <- E_7,4
res[14] <- E_7,5
res[15] <- E_7,6
res[16] <- E_7,7
res[17] <- E_7,8
res[18] <- E_7,9
res[19] <- E_7,10

Ordre 9
res[28] <- E_9,1
res[29] <- E_9,2
res[30] <- E_9,3
res[31] <- E_9,4
res[32] <- E_9,5
res[33] <- E_9,6
res[34] <- E_9,7
res[35] <- E_9,8
res[36] <- E_9,9

res[38] <- E_9,11
res[39] <- E_9,12
res[40] <- E_9,13
res[41] <- E_9,14
res[42] <- E_9,15
res[43] <- E_9,16
res[44] <- E_9,17
res[45] <- E_9,18
res[46] <- E_9,19
res[47] <- E_9,20
res[48] <- E_9,21
res[49] <- E_9,22
res[50] <- E_9,23
res[51] <- E_9,24
res[52] <- E_9,25
res[53] <- E_9,26
*/
#endif
