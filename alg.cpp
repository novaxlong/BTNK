//
// Created by 王恺忻 on 2018/11/12.
//

#include <queue>
#include <algorithm>
#include "alg.h"

using namespace std;

extern double **map;

Alg::Alg(VT &task, VW &worker) : taskList(task), workerList(worker), run(false) {
    bottleNeck = -INT_MAX;
    checked = vector<bool>(workerList.size(), false);
    workerMatchedList = VT(worker.size());
    taskMatchedList = VW(task.size());
}

void Alg::simpleGreedy() {
    VT taskWaiting;
    VW workerWaiting;
    double minValue;
    int minIndex;

    for (int i = 0, j = 0; i < taskList.size() || j < workerList.size(); ) {
        if (j == workerList.size()
                || (i != taskList.size() && taskList[i].nodeArr < workerList[j].nodeArr)) {
            if (!workerWaiting.empty()) {
                minValue = INT_MAX;
                minIndex = -1;
                for (int k = 0; k < workerWaiting.size(); ++k) {
                    if (map[taskList[i].nodeId][workerWaiting[k].nodeId] >=0
                            && map[taskList[i].nodeId][workerWaiting[k].nodeId] < minValue) {
                        minValue = map[taskList[i].nodeId][workerWaiting[k].nodeId];
                        minIndex = k;
                    }
                }
                if (minIndex != -1) {
                    taskMatchedList[taskList[i].nodeId] = workerWaiting[minIndex];
                    workerMatchedList[workerWaiting[minIndex].nodeId] = taskList[i];
                    bottleNeckPair = (bottleNeck > minValue ?
                                      bottleNeckPair : make_pair(taskList[i].nodeId, workerWaiting[minIndex].nodeId));
                    bottleNeck = (bottleNeck > minValue ? bottleNeck : minValue);
                    workerWaiting.erase(workerWaiting.begin() + minIndex);
                }
            } else
                taskWaiting.push_back(taskList[i]);
            i++;
        } else {
            if (!taskWaiting.empty()) {
                minValue = INT_MAX;
                minIndex = -1;
                for (int k = 0; k < taskWaiting.size(); ++k) {
                    if (map[taskWaiting[k].nodeId][workerList[j].nodeId] >=0
                            && map[taskWaiting[k].nodeId][workerList[j].nodeId] < minValue) {
                        minValue = map[taskWaiting[k].nodeId][workerList[j].nodeId];
                        minIndex = k;
                    }
                }
                if (minIndex != -1) {
                    taskMatchedList[taskWaiting[minIndex].nodeId] = workerList[j];
                    workerMatchedList[workerList[j].nodeId] = taskWaiting[minIndex];
                    bottleNeckPair = (bottleNeck > minValue ?
                            bottleNeckPair : make_pair(taskWaiting[minIndex].nodeId, workerList[j].nodeId));
                    bottleNeck = (bottleNeck > minValue ? bottleNeck : minValue);
                    taskWaiting.erase(taskWaiting.begin() + minIndex);
                }
            } else
                workerWaiting.push_back(workerList[j]);
            j++;
        }
    }
    outputResult("Greedy");
    run = true;
}

static int countdd = 0;

bool Alg::occupied(double d, task& t, VN& chain) {
    if (d < 0.37)
        countdd = countdd;
    countdd++;
//    printf("line 82----%d\n", countdd);
    if (map[t.nodeId][bottleNeckPair.second] < d && map[t.nodeId][bottleNeckPair.second] >= 0) {
        chain.push_back(workerList[bottleNeckPair.second]);
        chain.push_back(t);
        countdd--;
//        printf("line 87----%d\n", countdd);
        return true;
    }           // available

    for (auto& worker : workerList) {
        if (checked[worker.nodeId])
            continue;
        if (map[t.nodeId][worker.nodeId] < d && map[t.nodeId][worker.nodeId] >= 0) {
            checked[worker.nodeId] = true;
            task next = workerMatchedList[worker.nodeId];
            if (occupied(d, next, chain)) {
                chain.push_back(worker);
                chain.push_back(t);
                checked[worker.nodeId] = false;
                countdd--;
//                printf("line 102----%d\n", countdd);
                return true;
            }
            checked[worker.nodeId] = false;
        }
    }
    countdd--;
//    printf("line 109----%d\n", countdd);
    return false;
}

void Alg::swapChain() {
    if (!run) {
        printf("PLEASE run ANOTHER method to gain a matching.\n");
        return;
    }
    VN chain = VN();
    while (occupied(bottleNeck, taskList[bottleNeckPair.first], chain)) {
        for (int i = 0; i < chain.size(); ++i) {
            if (i & 1) {       // task
                taskMatchedList[chain[i].nodeId] = chain[i-1];
                workerMatchedList[chain[i-1].nodeId] = chain[i];
            } else {                // worker
                taskMatchedList[chain[i+1].nodeId] = chain[i];
                workerMatchedList[chain[i].nodeId] = chain[i+1];
            }
        }
        for (int i = 0; i < taskMatchedList.size(); ++i) {
            int workerId = taskMatchedList[i].nodeId;
            if (i == 0) {
                bottleNeckPair = make_pair(i, workerId);
                bottleNeck = map[i][workerId];
            } else {
                bottleNeckPair = (bottleNeck > map[i][workerId] ? bottleNeckPair : make_pair(i, workerId));
                bottleNeck = (bottleNeck > map[i][workerId] ? bottleNeck : map[i][workerId]);
            }
        }
        chain.clear();
    }
    outputResult("Optimal");
}

void Alg::qLearning() {
    outputResult("Q-Learning");
    // This is a new line of algorithm.
    run = true;
}

void Alg::outputResult(const char *algMethod) {
    printf("%s Performances: \n", algMethod);
    printf("bottleNeck = %lf\n", bottleNeck);
    printf("Tasks Matching\n");
    for (int i = 0; i < taskMatchedList.size(); ++i)
        printf("%4d", i);
    printf("\n");
    for (auto &worker : taskMatchedList)
        printf("%4d", worker.nodeId);
    printf("\nWorkers Matching\n");
    for (int i = 0; i < workerMatchedList.size(); ++i)
        printf("%4d", i);
    printf("\n");
    for (auto &task : workerMatchedList)
        printf("%4d", task.nodeId);
    printf("\nBottleNeckPair = <task(%d), worker(%d)>.\n\n", bottleNeckPair.first, bottleNeckPair.second);
}



