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
    DataProcess::splitSequence(nodeList, taskList, workerList);
    DataProcess::generateMatrix("../syntheticData/QLMatrix.dat", taskList, workerList);
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

void QL::Greedy(VT &taskList, VW &workerList, VW &taskMatchedList, VT &workerMatchedList) {
    VT taskWaiting;
    VW workerWaiting;
    double minValue;
    int minIndex;

    for (int i = 0, j = 0; i < taskList.size() || j < workerList.size(); ) {
        if (j == workerList.size()
            || (i != taskList.size() && taskList[i].nodeArr < workerList[j].nodeArr)) {
            // delete expired tasks
            for (int k = 0; k < taskWaiting.size(); ++k) {
                if (taskWaiting[k].nodeArr + taskWaiting[k].nodeDur < taskList[i].nodeArr) {
                    taskWaiting.erase(taskWaiting.begin() + k);
                    k--;
                }
            }
            // delete expired workers
            for (int k = 0; k < workerWaiting.size(); ++k) {
                if (workerWaiting[k].nodeArr + workerWaiting[k].nodeDur < taskList[i].nodeArr) {
                    workerWaiting.erase(workerWaiting.begin() + k);
                    k--;
                }
            }
            // if it exists workers waiting for a new task, these workers must have one to be assigned to the task
            if (!workerWaiting.empty()) {
                minValue = INT_MAX;
                minIndex = -1;
                for (int k = 0; k < workerWaiting.size(); ++k) {
                    if (QLMap[taskList[i].nodeId][workerWaiting[k].nodeId] >= 0 && QLMap[taskList[i].nodeId][workerWaiting[k].nodeId] < minValue) {
                        minValue = QLMap[taskList[i].nodeId][workerWaiting[k].nodeId];
                        minIndex = k;
                    }
                }
                taskMatchedList[taskList[i].nodeId] = workerWaiting[minIndex];
                workerMatchedList[workerWaiting[minIndex].nodeId] = taskList[i];
                bottleNeckPair = (bottleNeck > minValue ? bottleNeckPair : make_pair(taskList[i], workerWaiting[minIndex]));
                bottleNeck = (bottleNeck > minValue ? bottleNeck : minValue);
                workerWaiting.erase(workerWaiting.begin() + minIndex);
            }
            else
                taskWaiting.push_back(taskList[i]);
            i++;
        }
        else {
            // delete expired tasks
            for (int k = 0; k < taskWaiting.size(); ++k) {
                if (taskWaiting[k].nodeArr + taskWaiting[k].nodeDur < workerList[j].nodeArr) {
                    taskWaiting.erase(taskWaiting.begin() + k);
                    k--;
                }
            }
            // delete expired workers
            for (int k = 0; k < workerWaiting.size(); ++k) {
                if (workerWaiting[k].nodeArr + workerWaiting[k].nodeDur < workerList[j].nodeArr) {
                    workerWaiting.erase(workerWaiting.begin() + k);
                    k--;
                }
            }
            // if it exists tasks waiting for a new worker, these tasks must have one to be assigned to the worker
            if (!taskWaiting.empty()) {
                minValue = INT_MAX;
                minIndex = -1;
                for (int k = 0; k < taskWaiting.size(); ++k) {
                    if (QLMap[taskWaiting[k].nodeId][workerList[j].nodeId] >=0
                        && QLMap[taskWaiting[k].nodeId][workerList[j].nodeId] < minValue) {
                        minValue = QLMap[taskWaiting[k].nodeId][workerList[j].nodeId];
                        minIndex = k;
                    }
                }
                taskMatchedList[taskWaiting[minIndex].nodeId] = workerList[j];
                workerMatchedList[workerList[j].nodeId] = taskWaiting[minIndex];
                bottleNeckPair = (bottleNeck > minValue ? bottleNeckPair : make_pair(taskWaiting[minIndex], workerList[j]));
                bottleNeck = (bottleNeck > minValue ? bottleNeck : minValue);
                taskWaiting.erase(taskWaiting.begin() + minIndex);
            }
            else
                workerWaiting.push_back(workerList[j]);
            j++;
        }
    }
}

bool QL::BFS(double d, task bottleNeckTask, VW &taskMatchedList, VT &workerMatchedList, VN &chain) {
    QN nodeQueue = QN();
    bool updated = false;
    vector<bool> checkedWorker = vector<bool>(workerMatchedList.size(), false);
    vector<worker> taskParent = vector<worker>(taskMatchedList.size());
    vector<task> workerParent = vector<task>(workerMatchedList.size());

    nodeQueue.push_back(bottleNeckTask);
    while (!nodeQueue.empty()) {
        node firstItem = nodeQueue.front();
        if (firstItem.nodeType == 0) {                  // task
            for (node& worker : workerList) {
                if (!checkedWorker[worker.nodeId] && QLMap[firstItem.nodeId][worker.nodeId] < d && QLMap[firstItem.nodeId][worker.nodeId] >= 0) {
                    nodeQueue.push_back(worker);
                    workerParent[worker.nodeId] = firstItem;
                    checkedWorker[worker.nodeId] = true;
                }
            }
        }
        else {            // worker
            if (firstItem.nodeId == bottleNeckPair.second.nodeId) {
                updated = true;
                break;
            }
            else {
                nodeQueue.push_back(workerMatchedList[firstItem.nodeId]);
                taskParent[workerMatchedList[firstItem.nodeId].nodeId] = firstItem;
            }
        }
        nodeQueue.pop_front();
    }
    if (updated) {
        node lastTwoItem = bottleNeckPair.second;
        node lastItem = workerParent[lastTwoItem.nodeId];
        while (lastItem.nodeId != bottleNeckPair.first.nodeId) {
            chain.push_back(lastTwoItem);
            chain.push_back(lastItem);
            lastTwoItem = taskParent[lastItem.nodeId];
            lastItem = workerParent[lastTwoItem.nodeId];
        }
        chain.push_back(lastTwoItem);
        chain.push_back(lastItem);
        return true;
    }
    return false;
}

double QL::swapChain(VT &taskList, VW &workerList, VW &taskMatchedList, VT &workerMatchedList) {
    Greedy(taskList, workerList, taskMatchedList, workerMatchedList);
    VN chain = VN();
    while (BFS(bottleNeck, bottleNeckPair.first, taskMatchedList, workerMatchedList, chain)) {
        for (int i = 0; i < chain.size(); ++i) {
            if (i & 1) {       // task
                taskMatchedList[chain[i].nodeId] = chain[i-1];
                workerMatchedList[chain[i-1].nodeId] = chain[i];
            } else {                // worker
                taskMatchedList[chain[i+1].nodeId] = chain[i];
                workerMatchedList[chain[i].nodeId] = chain[i+1];
            }
        }
        bottleNeck = 0.0;
        for (node& task : taskList) {
            int workerId = taskMatchedList[task.nodeId].nodeId;
            bottleNeckPair = (bottleNeck > QLMap[task.nodeId][workerId] ? bottleNeckPair : make_pair(task, taskMatchedList[task.nodeId]));
            bottleNeck = (bottleNeck > QLMap[task.nodeId][workerId] ? bottleNeck : QLMap[task.nodeId][workerId]);
        }
        chain.clear();
    }
    return bottleNeck;
}

double QL::RQL(int len, int lambda, int upperBound) {
    init(len, lambda, upperBound);
    while ((int) taskList.size() > 0)
        taskList.erase(taskList.begin());
    while ((int) workerList.size() > 0)
        workerList.erase(workerList.begin());

    int count = 0, i = 0, size_l, size_r, last_l, last_r;
    double weight, btnk = 0.;
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
            VW taskMatchedList = VW((unsigned long) len);
            VT workerMatchedList = VT((unsigned long) len);
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
            Q[last_l][last_r][count - LOWER_BOUND][lt - LOWER_BOUND] += ALPHA * ((1.0 / weight) + temp - Q[last_l][last_r][count - LOWER_BOUND][lt - LOWER_BOUND]);
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
    finale(len);
    return btnk;
}


