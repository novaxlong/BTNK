//
// Created by 王恺忻 on 2018/11/10.
//

#ifndef BOTTLENECK_DATA_H
#define BOTTLENECK_DATA_H

#include <vector>

using namespace std;

struct node {
    int nodeType, nodeId, nodeArr, nodeDur;
    node() : nodeType(0), nodeId(0), nodeArr(0), nodeDur(0) {}
};

typedef node task;
typedef node worker;
typedef vector<task> VT;
typedef vector<worker> VW;
typedef vector<node> VN;

class DataProcess {

    /**
     * This class is going to generate and read different attributes of data from files.
     * generate(const char*, int) is to generate a data file in directory syntheticData/
     * and the data structure is a tuple, denoted by (type, id, arrive_time, duration)
     */

public:
    static void generateSequence(const char *fileName, int len, int lambda, int upperBound);
    static void readSequence(const char *fileName, VN& nodeList);
    static void splitSequence(VN nodeList, VT& taskList, VW& workerList);
    static void generateMatrix(const char *fileName, VT taskList, VW workerList);
    static void readMatrix(const char *fileName, int len, double* map[]);
    static void mallocMap(int len, double ***map);
    static void freeMap(int len, double ***map);
};

#endif //BOTTLENECK_DATA_H
