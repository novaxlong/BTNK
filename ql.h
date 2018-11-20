//
// Created by 王恺忻 on 2018/11/17.
//

#ifndef BOTTLENECK_QL_H
#define BOTTLENECK_QL_H

#include "data.h"
#include <list>

#define LOWER_BOUND 10
#define UPPER_BOUND 20
#define LEN (UPPER_BOUND - LOWER_BOUND + 1)
#define MAXS 50     // cannot initialize a big number. My computer does not have so much memory.
#define ALPHA 5e-2
#define INIT 1e-4

using namespace std;
typedef list<node> QN;

class QL {
public:
    pair<task, worker> bottleNeckPair;
    double bottleNeck;
    QL();
    double swapChain(VT& taskList, VW& workerList, VW& taskMatchedList, VT& workerMatchedList);
    int pick(int start, int duration);
    void init(int len, int lambda, int upperBound);
    void finale(int len);
    double RQL(int len, int lambda, int upperBound);

private:
    void Greedy(VT &taskList, VW &workerList, VW &taskMatchedList, VT &workerMatchedList);
    double **QLMap;
    VN nodeList;
    VT taskList;
    VW workerList;
    bool BFS(double d, task bottleNeckTask, VW &taskMatchedList, VT &workerMatchedList, VN& chain);
    double Q[MAXS][MAXS][LEN][LEN];
};


#endif //BOTTLENECK_QL_H
