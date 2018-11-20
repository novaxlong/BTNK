//
// Created by 王恺忻 on 2018/11/12.
//

#ifndef BOTTLENECK_ALG_H
#define BOTTLENECK_ALG_H

#include "data.h"
#include <list>

// parameters for QL
#define UPPER_BOUND 10
#define LOWER_BOUND 5
#define MAXSIZE 50
#define LEN (UPPER_BOUND - LOWER_BOUND + 1)
#define INIT 1e-4
#define ALPHA 5e-2

// parameters for CandGreedy
#define SIZE 10

typedef pair<int, int> PII;     // < L_index, R_index >

// Simple Greedy
double onlineGreedy(const VVD& cost, VVI& L, VVI& R, VI& Lmate, VI& Rmate, PII& btnkPair);

// Candidate Based Greedy
double CandGreedy(const VVD& cost, VVI& L, VVI& R, VI& Lmate, VI& Rmate, PII& btnkPair);

// Optimal
bool BFS(const VVD& cost, VVI& L, VVI& R, VI& Lmate, VI& Rmate, double btnk, PII& btnkPair, VI& chain);
double swapChain(const VVD& cost, VVI& L, VVI& R, VI& Lmate, VI& Rmate, PII& btnkPair);

// RQL
void initQ();
int tick(const VVI& seq, VVI& L, VVI& R, int start, int t);
void RQL(const VVD& cost, VVI& seq);
double QL(const VVD& cost, VVI& seq, VI& Lmate, VI& Rmate, PII& btnkPair);

#endif //BOTTLENECK_ALG_H
