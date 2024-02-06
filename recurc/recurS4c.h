/* 06-02-2024 */
/* alex */
/* recurS4c.h */
#ifndef _RECURS4C_H
#define _RECURS4C_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "esq.h"

void recS4_c(double complex * res, double complex x, int order);

/* Composició SS(S4)
Ordre 1
res[0] <- M_1,1

Ordre 5
res[1] <- M_5,1

Ordre 6
res[2] <- M_6,1

Ordre 7
res[3] <- M_7,1
res[4] <- M_7,2

Ordre 8
res[5] <- M_8,1
res[6] <- M_8,2

Ordre 9
res[7] <- M_9,1
res[8] <- M_9,2
res[9] <- M_9,3

Ordre 10
res[10] <- M_10,1
res[11] <- M_10,2
res[12] <- M_10,3

Ordre 11
res[13] <- M_11,1
res[14] <- M_11,2
res[15] <- M_11,3
res[16] <- M_11,4
res[17] <- M_11,5

Ordre 12
res[18] <- M_11,1
res[19] <- M_11,2
res[20] <- M_11,3
res[21] <- M_11,4
res[22] <- M_11,5
res[23] <- M_11,6

Ordre 13
res[24] <- M_11,1
res[25] <- M_11,2
res[26] <- M_11,3
res[27] <- M_11,4
res[28] <- M_11,5
res[29] <- M_11,6
res[30] <- M_11,7
res[31] <- M_11,8
res[32] <- M_11,9

Ordre 14
res[33] <- M_11,1
res[34] <- M_11,2
res[35] <- M_11,3
res[36] <- M_11,4
res[37] <- M_11,5
res[38] <- M_11,6
res[39] <- M_11,7
res[40] <- M_11,8
res[41] <- M_11,9
res[42] <- M_11,10
res[43] <- M_11,11

Ordre 15
res[44] <- M_11,1
res[45] <- M_11,2
res[46] <- M_11,3
res[47] <- M_11,4
res[48] <- M_11,5
res[49] <- M_11,6
res[50] <- M_11,7
res[51] <- M_11,8
res[52] <- M_11,9
res[53] <- M_11,10
res[54] <- M_11,11
res[55] <- M_11,12
res[56] <- M_11,13
res[57] <- M_11,14
res[58] <- M_11,15
res[59] <- M_11,16

*/
#endif
