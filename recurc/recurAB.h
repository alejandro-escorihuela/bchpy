/* 28-02-2022 */
/* alex */
/* recurAB.h */
#ifndef _RECURAB_H
#define _RECURAB_H
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "esq.h"

void recAsim(double * res, double x, int order, int rkn);
void recBsim(double * res, double x, int order, int rkn);
void recA(double * res, double x, int order, int rkn);
void recB(double * res, double x, int order, int rkn);

/*
  Escissió simètrica
  Ordre 1             Ordre 9 
  res[0] <- E_1,1     res[28] <- E_9,1    res[59] <- E_9,32
  res[1] <- E_1,2     res[29] <- E_9,2    res[60] <- E_9,33
  Ordre 3             res[30] <- E_9,3    res[61] <- E_9,34
  res[2] <- E_3,1     res[31] <- E_9,4    res[62] <- E_9,35
  res[3] <- E_3,2     res[32] <- E_9,5    res[63] <- E_9,36
  Ordre 5             res[33] <- E_9,6    res[64] <- E_9,37
  res[4] <- E_5,1     res[34] <- E_9,7    res[65] <- E_9,38
  res[5] <- E_5,2     res[35] <- E_9,8    res[66] <- E_9,39
  res[6] <- E_5,3     res[36] <- E_9,9    res[67] <- E_9,40
  res[7] <- E_5,4     res[37] <- E_9,10   res[68] <- E_9,41
  res[8] <- E_5,5     res[38] <- E_9,11   res[69] <- E_9,42
  res[9] <- E_5,6     res[39] <- E_9,12   res[70] <- E_9,43
  Ordre 7             res[40] <- E_9,13   res[71] <- E_9,44
  res[10] <- E_7,1    res[41] <- E_9,14   res[72] <- E_9,45
  res[11] <- E_7,2    res[42] <- E_9,15   res[73] <- E_9,46
  res[12] <- E_7,3    res[43] <- E_9,16   res[74] <- E_9,47
  res[13] <- E_7,4    res[44] <- E_9,17   res[75] <- E_9,48
  res[14] <- E_7,5    res[45] <- E_9,18   res[76] <- E_9,49
  res[15] <- E_7,6    res[46] <- E_9,19   res[77] <- E_9,50
  res[16] <- E_7,7    res[47] <- E_9,20   res[78] <- E_9,51
  res[17] <- E_7,8    res[48] <- E_9,21   res[79] <- E_9,52
  res[18] <- E_7,9    res[49] <- E_9,22   res[80] <- E_9,53
  res[19] <- E_7,10   res[50] <- E_9,23   res[81] <- E_9,54
  res[20] <- E_7,11   res[51] <- E_9,24   res[82] <- E_9,55
  res[21] <- E_7,12   res[52] <- E_9,25   res[83] <- E_9,56
  res[22] <- E_7,13   res[53] <- E_9,26 
  res[23] <- E_7,14   res[54] <- E_9,27 
  res[24] <- E_7,15   res[55] <- E_9,28 
  res[25] <- E_7,16   res[56] <- E_9,29 
  res[26] <- E_7,17   res[57] <- E_9,30 
  res[27] <- E_7,18   res[58] <- E_9,31 

*/
/*
  Escissió no simètrica
  Ordre 1            res[31] <- E_7,9    res[68] <- E_8,28   res[105] <- E_9,35
  res[0] <- E_1,1    res[32] <- E_7,10   res[69] <- E_8,29   res[106] <- E_9,36
  res[1] <- E_1,2    res[33] <- E_7,11   res[70] <- E_8,30   res[107] <- E_9,37
  Ordre 2            res[34] <- E_7,12   Ordre 9             res[108] <- E_9,38
  res[2] <- E_2,1    res[35] <- E_7,13   res[71] <- E_9,1    res[109] <- E_9,39
  Ordre 3            res[36] <- E_7,14   res[72] <- E_9,2    res[110] <- E_9,40
  res[3] <- E_3,1    res[37] <- E_7,15   res[73] <- E_9,3    res[111] <- E_9,41
  res[4] <- E_3,2    res[38] <- E_7,16   res[74] <- E_9,4    res[112] <- E_9,42
  Ordre 4            res[39] <- E_7,17   res[75] <- E_9,5    res[113] <- E_9,43
  res[5] <- E_4,1    res[40] <- E_7,18   res[76] <- E_9,6    res[114] <- E_9,44
  res[6] <- E_4,2    Ordre 8             res[77] <- E_9,7    res[115] <- E_9,45
  res[7] <- E_4,3    res[41] <- E_8,1    res[78] <- E_9,8    res[116] <- E_9,46
  Ordre 5            res[42] <- E_8,2    res[79] <- E_9,9    res[117] <- E_9,47
  res[8] <- E_5,1    res[43] <- E_8,3    res[80] <- E_9,10   res[118] <- E_9,48
  res[9] <- E_5,2    res[44] <- E_8,4    res[81] <- E_9,11   res[119] <- E_9,49
  res[10] <- E_5,3   res[45] <- E_8,5    res[82] <- E_9,12   res[120] <- E_9,50
  res[11] <- E_5,4   res[46] <- E_8,6    res[83] <- E_9,13   res[121] <- E_9,51
  res[12] <- E_5,5   res[47] <- E_8,7    res[84] <- E_9,14   res[122] <- E_9,52
  res[13] <- E_5,6   res[48] <- E_8,8    res[85] <- E_9,15   res[123] <- E_9,53
  Ordre 6            res[49] <- E_8,9    res[86] <- E_9,16   res[124] <- E_9,54
  res[14] <- E_6,1   res[50] <- E_8,10   res[87] <- E_9,17   res[125] <- E_9,55
  res[15] <- E_6,2   res[51] <- E_8,11   res[88] <- E_9,18   res[126] <- E_9,56
  res[16] <- E_6,3   res[52] <- E_8,12   res[89] <- E_9,19 
  res[17] <- E_6,4   res[53] <- E_8,13   res[90] <- E_9,20 
  res[18] <- E_6,5   res[54] <- E_8,14   res[91] <- E_9,21 
  res[19] <- E_6,6   res[55] <- E_8,15   res[92] <- E_9,22 
  res[20] <- E_6,7   res[56] <- E_8,16   res[93] <- E_9,23 
  res[21] <- E_6,8   res[57] <- E_8,17   res[94] <- E_9,24 
  res[22] <- E_6,9   res[58] <- E_8,18   res[95] <- E_9,25 
  Ordre 7            res[59] <- E_8,19   res[96] <- E_9,26 
  res[23] <- E_7,1   res[60] <- E_8,20   res[97] <- E_9,27 
  res[24] <- E_7,2   res[61] <- E_8,21   res[98] <- E_9,28 
  res[25] <- E_7,3   res[62] <- E_8,22   res[99] <- E_9,29 
  res[26] <- E_7,4   res[63] <- E_8,23   res[100] <- E_9,30 
  res[27] <- E_7,5   res[64] <- E_8,24   res[101] <- E_9,31 
  res[28] <- E_7,6   res[65] <- E_8,25   res[102] <- E_9,32 
  res[29] <- E_7,7   res[66] <- E_8,26   res[103] <- E_9,33 
  res[30] <- E_7,8   res[67] <- E_8,27   res[104] <- E_9,34

*/


#endif
