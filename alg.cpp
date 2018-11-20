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
    checkedWorker = vector<bool>(workerList.size(), false);
    checkedTask = vector<bool>(taskList.size(), false);
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
    for (int l = 0; l < workerWaiting.size(); ++l) {
        printf("Discard worker[%d].\n", workerWaiting[l].nodeId);
    }
    for (int m = 0; m < taskWaiting.size(); ++m) {
        printf("Discarding task[%d].\n", taskWaiting[m].nodeId);
    }
    outputResult("Greedy");
    run = true;
}

bool Alg::IBFS(double d, task bottleNeckTask, VN& chain) {
    QN nodeQueue = QN();
    bool updated = false;
    int *taskParent, *workerParent;
    taskParent = (int*) malloc(taskList.size() * sizeof(int));
    workerParent = (int*) malloc(workerList.size() * sizeof(int));

    nodeQueue.push_back(bottleNeckTask);
    while (!nodeQueue.empty()) {
//        printf("%ld\n", nodeQueue.size());
        node firstItem = nodeQueue.front();
        if (firstItem.nodeType == 0) {                  // task
            for (node &worker : workerList) {
                if (!checkedWorker[worker.nodeId] && map[firstItem.nodeId][worker.nodeId] < d && map[firstItem.nodeId][worker.nodeId] >= 0) {
                    nodeQueue.push_back(worker);
                    workerParent[worker.nodeId] = firstItem.nodeId;
                    checkedWorker[worker.nodeId] = true;
                }
            }
        } else {            // worker
            if (firstItem.nodeId == bottleNeckPair.second) {
                updated = true;
                break;
            }
            else {
//                if (!checkedTask[workerMatchedList[firstItem.nodeId].nodeId]) {
//                    checkedTask[workerMatchedList[firstItem.nodeId].nodeId] = true;
                    nodeQueue.push_back(workerMatchedList[firstItem.nodeId]);
                    taskParent[workerMatchedList[firstItem.nodeId].nodeId] = firstItem.nodeId;
//                }

            }
        }
        nodeQueue.pop_front();
    }
    if (updated) {
        node lastTwoItem = workerList[bottleNeckPair.second];
        node lastItem = taskList[workerParent[bottleNeckPair.second]];
        while (lastItem.nodeId != bottleNeckPair.first) {
            chain.push_back(lastTwoItem);
            chain.push_back(lastItem);
            lastTwoItem = workerList[taskParent[lastItem.nodeId]];
            lastItem = taskList[workerParent[lastTwoItem.nodeId]];
        }
        chain.push_back(lastTwoItem);
        chain.push_back(lastItem);
        free(taskParent);
        free(workerParent);
        return true;
    }
    free(taskParent);
    free(workerParent);
    return false;
}

void Alg::swapChain() {
    if (!run) {
        printf("PLEASE run ANOTHER method to gain a matching.\n");
        return;
    }
    VN chain = VN();
    while (IBFS(bottleNeck, taskList[bottleNeckPair.first], chain)) {
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
        checkedWorker = vector<bool>(workerList.size(), false);
        checkedTask = vector<bool>(taskList.size(), false);
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



