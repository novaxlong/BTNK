#include <iostream>
#include <vector>
#include <cstdlib>
#include "data.h"
#include "alg.h"
//#include "ql.h"

using namespace std;

int main(int argc, const char* argv[]) {

    int len = 20, lambda = 3, lowerBound = 10, upperBound = 20;
    double btnk;
    VVI seq, L, R;
    VVD cost;

    generateSequence("../syntheticData/sequence.dat", len, lambda, lowerBound, upperBound);
    readSequence("../syntheticData/sequence.dat", seq);
    splitSequence(seq, L, R);
    generateMatrix("../syntheticData/matrix.dat", L, R);
    readMatrix("../syntheticData/matrix.dat", len, cost);

    VI Lmate, Rmate;
    PII pair;

    VVI LL, RR;
    for (int i = 0; i < 11; ++i) {
        LL.push_back(L[i+5]);
    }
    for (int i = 0; i < 7; ++i) {
        RR.push_back(R[i+5]);
    }

    btnk = onlineGreedy(cost, LL, RR, Lmate, Rmate, pair);
    printf("%.3f: < %d, %d >\n", btnk, LL[pair.first][1], RR[pair.second][1]);
    printf("\n\n");

    btnk = swapChain(cost, LL, RR, Lmate, Rmate, pair);
    printf("%.3f: < %d, %d >\n", btnk, LL[pair.first][1], RR[pair.second][1]);
    printf("\n\n");
//
//    Alg algorithm(taskList, workerList);
//
//    algorithm.simpleGreedy();
//    algorithm.swapChain();
//
//    QL myQL;
//    myQL.RQL(len, lambda, upperBound);

    return 0;
}

























































































//using namespace std;
//
//int taskNum, workerNum;
//int map[10000][10000];
//bool used[10000] = { false };
//bool checked[10000] = { false };
//vector<task> tasks;
//vector<worker> workers;
//
//void generate(int n, const char* fileName) {
//    FILE *fp = fopen(fileName, "wb");
//    if (fp == nullptr) {
//        cout << "File does not exists." << endl;
//        exit(-1);
//    }
//    fprintf(fp, "%d\n", n);
//    srand((unsigned int) time(nullptr));
//    for (int i = 0; i < n; ++i) {
//        for (int j = 0; j < n; ++j) {
//            fprintf(fp, "%d %d %d\n", i, j, rand() % 20);
//        }
//    }
//    fclose(fp);
//}
//
//void readData(const char* fileName) {
//    /*
//     * The function read data from bottleneck.txt.
//     * The first line of the bottleneck.txt is a number.
//     * The remaining lines of the rest of this file contain the information of distances
//     * between different pairs. We use this function to restore.
//     */
//    FILE *fp = fopen(fileName, "rb");
//    if (fp == nullptr) {
//        cout << "File does not exists." << endl;
//        exit(-1);
//    }
//    int n, i, j, distance;
//    fscanf(fp, "%d", &n);
//    taskNum = workerNum = n;
//    while (fscanf(fp, "%d %d %d", &i, &j, &distance) == 3) {
//        map[i][j] = distance;
//    }
//    fclose(fp);
//}
//
//void initMatch() {
//    tasks = vector<task>((unsigned long) taskNum);
//    workers = vector<worker>((unsigned long) workerNum);
//    for (int i = 0; i < taskNum; ++i) {
//        tasks[i].index = workers[i].index = i;
//        tasks[i].matchedIndex = workers[i].matchedIndex = i;
//        tasks[i].distance = workers[i].distance = map[i][i];
//        used[workers[i].index] = true;
//    }
//}
//
//bool available(int d, task t, vector<node>& chain) {
//    for (int i = 0; i < workers.size(); i++) {
//        if (checked[workers[i].index] || used[workers[i].index])
//            continue;
//        if (map[t.index][workers[i].index] < d) {
//            chain.push_back(workers[i]);
//            chain.push_back(t);
//            return true;
//        }
//    }
//    return false;
//}
//
//bool occupied(int d, task t, vector<node>& chain) {
//    for (int i = 0; i < workers.size(); ++i) {
//        if (checked[workers[i].index])
//            continue;
//        if (map[t.index][workers[i].index] < d) {
//            checked[workers[i].index] = true;
//            task nextTask = tasks[workers[i].matchedIndex];
//            if (available(d, nextTask, chain) || occupied(d, nextTask, chain)) {
//                chain.push_back(workers[i]);
//                chain.push_back(t);
//                checked[workers[i].index] = false;
//                return true;
//            }
//            checked[workers[i].index] = false;
//        }
//    }
//
//    return false;
//}
//
//int getBottleNeckValue(int& taskIndex, int& workerIndex) {
//    int bottleNeck = -1;
//    for (int i = 0; i < tasks.size(); ++i) {
//        if (bottleNeck < tasks[i].distance) {
//            bottleNeck = tasks[i].distance;
//            taskIndex = tasks[i].index;
//            workerIndex = tasks[i].matchedIndex;
//        }
//    }
//    return bottleNeck;
//}
//
//void output(int bottleNeck) {
//    cout << "Bottle Neck Value = " << bottleNeck << endl;
//    for (int i = 0; i < tasks.size(); ++i) {
//        cout << "task[" << i << "] matches worker[" << workers[tasks[i].matchedIndex].index << "]." << endl;
//    }
//    cout << endl;
//}
//
//void swapChain() {
//    int bottleNeckValue = -1, bottleNeckValueNew;
//    int bottleNeckTaskIndex, bottleNeckWorkerIndex;
//    while ((bottleNeckValueNew = getBottleNeckValue(bottleNeckTaskIndex, bottleNeckWorkerIndex)) != bottleNeckValue) {
//        output(bottleNeckValueNew);
//        bottleNeckValue = bottleNeckValueNew;
//        used[bottleNeckWorkerIndex] = false;
//        vector<node> chain = vector<node>();
//        occupied(bottleNeckValue, tasks[bottleNeckTaskIndex], chain);
//        for (int i = 0; i < chain.size(); ++i) {
//            if (i % 2 != 0) {
//                tasks[chain[i].index].matchedIndex = workers[chain[i-1].index].index;
//                tasks[chain[i].index].distance = map[tasks[chain[i].index].index][chain[i-1].index];
//            } else {
//                workers[chain[i].index].matchedIndex = tasks[chain[i+1].index].index;
//                workers[chain[i].index].distance = map[tasks[chain[i+1].index].index][chain[i].index];
//            }
//        }
//        used[bottleNeckWorkerIndex] = true;
//    }
//}
//
//int main(int argc, const char* argv[]) {
//    generate(5, "../bottleneck.txt");
//    readData("../bottleneck.txt");
//    initMatch();
//    swapChain();
//
//    count();
//}