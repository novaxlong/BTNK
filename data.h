//
// Created by 王恺忻 on 2018/11/10.
//

#ifndef BOTTLENECK_DATA_H
#define BOTTLENECK_DATA_H

#include <vector>

using namespace std;

typedef vector<int> VI;
typedef vector<VI> VVI;
typedef vector<double> VD;
typedef vector<VD> VVD;

void generateSequence(const char *fileName, int len, int lambda, int lowerBound, int upperBound);
void readSequence(const char *fileName, VVI& seq);
void splitSequence(VVI& seq, VVI& L, VVI& R);
void generateMatrix(const char *fileName, VVI& L, VVI& R);
void readMatrix(const char *fileName, int len, VVD& mat);

#endif //BOTTLENECK_DATA_H
