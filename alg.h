//
// Created by 王恺忻 on 2018/11/12.
//

#ifndef BOTTLENECK_ALG_H
#define BOTTLENECK_ALG_H

#include "data.h"

class Alg {

    /**
     * This class is going to define some algorithms to solve the bottleneck matching
     * problem, including Greedy, SwapChain, Q_Learning. Besides, I will introduce some
     * new algorithms to this class if I come up with new ideas. Let's begin!
     */

public:
    double bottleNeck;
    pair<int, int> bottleNeckPair;      // < task_id, worker_id >
    VT taskList, workerMatchedList;
    VW workerList, taskMatchedList;

    Alg (VT& task, VW& worker);
    void simpleGreedy();              // use Greedy algorithm to create the matching
    void swapChain();           // optimal algorithm
    void qLearning();           // inspired by Mr. Wang

private:
    bool run;
    vector<bool> checked;
    void outputResult(const char *algMethod);
    bool occupied(double d, task& t, VN& chain);
};


#endif //BOTTLENECK_ALG_H
