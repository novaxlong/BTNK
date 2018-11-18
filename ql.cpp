//
// Created by 王恺忻 on 2018/11/17.
//

#include "ql.h"
#include <cstdlib>

QL::QL() {
    for (int i = 0; i < MAXS; ++i) {
        for (int j = 0; j < MAXS; ++j) {
            for (int k = 0; k < LEN; ++k) {
                for (int l = k; l < LEN; ++l) {
                    Q[i][j][k][l] = INIT * (double(rand()) / double(RAND_MAX));
                }
            }
        }
    }
    nodeList = VN();
    taskList = VT();
    workerList = VW();
}

void QL::init(int len, int lambda, int upperBound) {
    // Create the "training data" for Q-value
    DataProcess::generateSequence("../syntheticData/QLSequence.dat", len, lambda, upperBound);
    DataProcess::readSequence("../syntheticData/QLSequence.dat", nodeList);
    DataProcess::mallocMap(len, &QLMap);
    DataProcess::readMatrix("../syntheticData/QLMatrix.dat", len, QLMap);
}

void QL::finale(int len) {
    // Free the allocated memory by the function init.
    DataProcess::freeMap(len, &QLMap);
}

int QL::pick(int start, int duration) {
    int time = nodeList[start].nodeArr, i = start;
    int end_time = time + duration;
    while ((time = nodeList[i].nodeArr) < end_time) {
        if (nodeList[i].nodeType == 0)
            taskList.push_back(nodeList[i]);
        else
            workerList.push_back(nodeList[i]);
        if ((int) taskList.size() >= MAXS)
            taskList.erase(taskList.begin());
        if ((int) workerList.size() >= MAXS)
            workerList.erase(workerList.begin());
        int taskSize = (int) taskList.size();
        int workerSize = (int) workerList.size();
        for (int j = 0; j < taskSize; ++j) {
            if (time > taskList[j].nodeArr + taskList[j].nodeDur || taskList[j].nodeType == -1) {
                taskList.erase(taskList.begin() + j);
                j--;
                taskSize--;
            }
        }
        for (int j = 0; j < workerSize; ++j) {
            if (time > workerList[j].nodeArr + workerList[j].nodeDur || workerList[j].nodeType == -1) {
                workerList.erase(workerList.begin() + j);
                j--;
                workerSize--;
            }
        }
        i++;
        if (i >= (int) nodeList.size())
            return i;
    }
    return i;
}

double QL::swapChain(VT &taskList, VW &workerList, VW &taskMatchedList, VT &workerMatchedList) {
    return 0.;
}

double QL::RQL(int len, int lambda, int upperBound) {
    while ((int) taskList.size() > 0)
        taskList.erase(taskList.begin());
    while ((int) workerList.size() > 0)
        workerList.erase(workerList.begin());

    int count = 0, i = 0, size_l, size_r, last_l, last_r;
    double weight, btnk = 1.;
    i = pick(0, LOWER_BOUND);
    count = LOWER_BOUND;

    while (i < nodeList.size()) {
        int time = nodeList[i].nodeArr;
        size_l = (int) taskList.size();
        size_r = (int) workerList.size();
        int lt = count;
        double temp = Q[size_l][size_r][count - LOWER_BOUND][count - LOWER_BOUND];
        for (int j = count; j <= LOWER_BOUND; ++j) {
            if (Q[size_l][size_r][count - LOWER_BOUND][count - LOWER_BOUND] > temp) {
                temp = Q[size_l][size_r][count - LOWER_BOUND][count - LOWER_BOUND];
                lt = j;
            }
        }
        if (lt == count) {
            VW taskMatchedList = VW(workerList.size());
            VT workerMatchedList = VT(taskList.size());
            weight = swapChain(taskList, workerList, taskMatchedList, workerMatchedList);

            for (int j = 0; j < (int) taskMatchedList.size(); ++j) {
                if (taskMatchedList[i].nodeId == -1)
                    continue;
                else
                    taskList[j].nodeType = -1;
            }
            for (int j = 0; j < (int) workerMatchedList.size(); ++j) {
                if (workerMatchedList[i].nodeType == -1)
                    continue;
                else
                    workerList[j].nodeType = -1;
            }
            btnk = (weight < btnk) ? btnk : weight;
            size_l = (int) taskList.size();
            size_r = (int) workerList.size();
            for (int j = 0; j < size_l; ++j) {
                if (time > taskList[j].nodeArr + taskList[j].nodeDur || taskList[j].nodeType == -1) {
                    taskList.erase(taskList.begin() + j);
                    j--;
                    size_l--;
                }
            }
            for (int j = 0; j < size_r; ++j) {
                if (time > workerList[j].nodeArr + workerList[j].nodeDur || workerList[j].nodeType == -1) {
                    workerList.erase(workerList.begin() + j);
                    j--;
                    size_r--;
                }
            }
            last_l = size_l;
            last_r = size_r;
            i = pick(i, LOWER_BOUND);
            size_l = (int) taskList.size();
            size_r = (int) workerList.size();
            temp = 0.;
            for (int j = 0; j < LEN; ++j) {
                if (Q[size_l][size_r][0][j] > temp)
                    temp = Q[size_l][size_r][0][j];
            }
            Q[last_l][last_r][count - LOWER_BOUND][lt - LOWER_BOUND] += ALPHA * ((1.0 / btnk) + temp - Q[last_l][last_r][count - LOWER_BOUND][lt - LOWER_BOUND]);
            count = LOWER_BOUND;
        } else {
            last_l = size_l;
            last_r = size_r;
            i = pick(i, 1);
            size_l = (int) taskList.size();
            size_r = (int) workerList.size();
            temp = 0.;
            for (int j = 0; j < LEN; ++j) {
                if (Q[size_l][size_r][count - LOWER_BOUND + 1][j] > temp)
                    temp = Q[size_l][size_r][count - LOWER_BOUND + 1][j];
            }
            Q[last_l][last_r][count - LOWER_BOUND][lt - LOWER_BOUND] += ALPHA * (temp - Q[last_l][last_r][count - LOWER_BOUND][lt - LOWER_BOUND]);
            count += 1;
        }
    }
    return btnk;
}


