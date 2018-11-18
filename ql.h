//
// Created by 王恺忻 on 2018/11/17.
//

#ifndef BOTTLENECK_QL_H
#define BOTTLENECK_QL_H

#include "data.h"

#define LOWER_BOUND 10
#define UPPER_BOUND 20
#define LEN (UPPER_BOUND - LOWER_BOUND + 1)
#define MAXS 100
#define ALPHA 5e-2
#define INIT 1e-4


class QL {
public:
    QL();
    double swapChain(VT& taskList, VW& workerList, VW& taskMatchedList, VT& workerMatchedList);
    int pick(int start, int duration);
    void init(int len, int lambda, int upperBound);
    void finale(int len);
    double RQL(int len, int lambda, int upperBound);

private:
    double **QLMap;
    VN nodeList;
    VT taskList;
    VW workerList;
    double Q[MAXS][MAXS][LEN][LEN];
};


#endif //BOTTLENECK_QL_H
