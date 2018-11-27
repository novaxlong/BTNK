//
// Created by 王恺忻 on 2018/11/10.
//

#ifndef BOTTLENECK_DATA_H
#define BOTTLENECK_DATA_H

#include <vector>

/*
 * Modified on 2018.11.22 by Wang
 * We add 4 hyper parameters to generate the sequence interval variables and distance variables.
 * We modified the FUNCTION "generateSequence", "generateMatrix", "readMatrix".
 */
#define LOWER_BOUND_SEQ 10
#define UPPER_BOUND_SEQ 20
#define LOWER_BOUND_DIS 10
#define UPPER_BOUND_DIS 20

using namespace std;

typedef vector<int> VI;
typedef vector<VI> VVI;
typedef vector<double> VD;
typedef vector<VD> VVD;
typedef vector<bool> VB;

void generateSequence(const char *fileName, int len, int lambda);
void readSequence(const char *fileName, VVI& seq);
void splitSequence(VVI& seq, VVI& L, VVI& R);
void generateMatrix(const char *fileName, VVI& L, VVI& R);
void readMatrix(const char *fileName, int len, VVD& mat);

#endif //BOTTLENECK_DATA_H
