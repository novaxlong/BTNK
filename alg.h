//
// Created by 王恺忻 on 2018/11/12.
//

#ifndef BOTTLENECK_ALG_H
#define BOTTLENECK_ALG_H

#include "data.h"
#include <list>

typedef pair<int, int> PII;     // < L_index, R_index >

double onlineGreedy(const VVD& cost, VVI& L, VVI& R, VI& Lmate, VI& Rmate, PII& btnkPair);
double swapChain(const VVD& cost, VVI& L, VVI& R, VI& Lmate, VI& Rmate, PII& btnkPair);

#endif //BOTTLENECK_ALG_H
