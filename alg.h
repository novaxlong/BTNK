//
// Created by 王恺忻 on 2018/11/12.
//

#ifndef BOTTLENECK_ALG_H
#define BOTTLENECK_ALG_H

#include "data.h"
#include <list>


#define UPPER_BOUND 20
#define LOWER_BOUND 10
#define MAXSIZE 50
#define LEN (UPPER_BOUND - LOWER_BOUND + 1)
#define INIT 1e-4
#define ALPHA 5e-2

typedef pair<int, int> PII;     // < L_index, R_index >


double onlineGreedy(const VVD& cost, VVI& L, VVI& R, VI& Lmate, VI& Rmate, PII& btnkPair);
double swapChain(const VVD& cost, VVI& L, VVI& R, VI& Lmate, VI& Rmate, PII& btnkPair);

void initQ();
int tick(const VVI& seq, VVI& L, VVI& R, int start, int t);
double RQL(const VVD& cost, VVI& seq, VI& Lmate, VI& Rmate);

#endif //BOTTLENECK_ALG_H
