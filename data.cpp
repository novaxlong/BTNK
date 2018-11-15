//
// Created by 王恺忻 on 2018/11/10.
//

#include <random>
#include "data.h"

using namespace std;

void DataProcess::generateSequence(const char *fileName, int len, int lambda, int upperBound) {
    FILE *fp = fopen(fileName, "w");
    if (fp == nullptr) {
        printf("Cannot open file %s\n.", fileName);
        exit(-1);
    }
    if (len <= 0) {
        printf("Invalid length parameter %d\n.", len);
        exit(-2);
    }

    int num = 0, totalLen = 2 * len, taskLen = 0, workerLen = 0;
    int taskNodeId = 0, workerNodeId = 0, nodeArr = 0, nodeDur;
    random_device rd;
    mt19937 gen(rd());
    poisson_distribution<> pd(lambda);
    uniform_int_distribution<> uid(1, upperBound);

    while (num < totalLen) {
        int k = pd(gen), task_k, worker_k;
        if (num + k > totalLen)
            k = totalLen - num;
        if (taskLen >= workerLen) {
            task_k = (int) (k * ((double) (len - taskLen) / (totalLen - taskLen - workerLen)));
            worker_k = k - task_k;
        }
        else {
            worker_k = (int) (k * ((double) (len - workerLen) / (totalLen - taskLen - workerLen)));
            task_k = k - worker_k;
        }
        for (int i = 0; i < task_k; ++i) {
            nodeDur = uid(gen);
            fprintf(fp, "0 %d %d %d\n", taskNodeId++, nodeArr, nodeDur);
        }
        for (int j = 0; j < worker_k; ++j) {
            nodeDur = uid(gen);
            fprintf(fp, "1 %d %d %d\n", workerNodeId++, nodeArr, nodeDur);
        }
        taskLen += task_k;
        workerLen += worker_k;
        nodeArr++;
        num += k;
    }
    fclose(fp);
}

void DataProcess::readSequence(const char *fileName, vector<node>& nodeList) {
    FILE *fp = fopen(fileName, "r");
    if (fp == nullptr) {
        printf("Cannot open file %s\n.", fileName);
        exit(-1);
    }

    int nodeType, nodeId, nodeArr, nodeDur;
    node temp;
    while (fscanf(fp, "%d %d %d %d", &nodeType, &nodeId, &nodeArr, &nodeDur) == 4) {
        temp.nodeType = nodeType;
        temp.nodeId = nodeId;
        temp.nodeArr = nodeArr;
        temp.nodeDur = nodeDur;
        nodeList.push_back(temp);
    }
    fclose(fp);
}

void DataProcess::splitSequence(VN nodeList, VT &taskList, VW& workerList) {
    taskList.clear();
    workerList.clear();
    for (node &curNode : nodeList) {
        if (curNode.nodeType == 0)
            taskList.push_back(curNode);
        else
            workerList.push_back(curNode);
    }
}

void DataProcess::generateMatrix(const char *fileName, VT taskList, VW workerList) {
    FILE *fp = fopen(fileName, "w");
    if (fp == nullptr) {
        printf("Cannot open file %s\n.", fileName);
        exit(-1);
    }

    double bottleNeck;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> urd(0, 1);

    for (auto &curTask : taskList) {
        for (auto &curWorker : workerList) {
            if ((curTask.nodeArr - (curWorker.nodeArr + curWorker.nodeDur))
                    * ((curTask.nodeArr + curTask.nodeDur) - curWorker.nodeArr) <= 0) {
                bottleNeck = urd(gen);
                fprintf(fp, "%lf ", bottleNeck);
            }
            else {
                bottleNeck = -1;
                fprintf(fp, "%d ", (int) bottleNeck);
            }
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

void DataProcess::readMatrix(const char *fileName, int len, double **map) {
    FILE *fp = fopen(fileName, "r");
    if (fp == nullptr) {
        printf("Cannot open file %s\n.", fileName);
        exit(-1);
    }
    if (map == nullptr) {
        map = (double **) malloc(len * sizeof(double *));
        for (int i = 0; i < len; ++i)
            map[i] = (double *) malloc(len * sizeof(double));
    }

    double value;
    for (int i = 0; i < len; ++i) {
        for (int j = 0; j < len; ++j) {
            fscanf(fp, "%lf", &value);
            map[i][j] = value;
        }
    }

    fclose(fp);
}

