#include <iostream>
#include <vector>
#include <cstdlib>
#include "data.h"
#include "alg.h"

using namespace std;



int main(int argc, const char* argv[]) {

    int len = 100, lambda = 3;

//    train function Q.
    initQ();
//    for (int i = 0; i < 100000; ++i) {
//        int QLLen = len, QLLambda = lambda;
//        VVI QLSeq, QLL, QLR;
//        VVD QLCost;
//        VI QLLmate, QLRmate;
//        generateSequence("../syntheticData/QLSequence.dat", QLLen, QLLambda);
//        readSequence("../syntheticData/QLSequence.dat", QLSeq);
//        splitSequence(QLSeq, QLL, QLR);
//        generateMatrix("../syntheticData/QLMatrix.dat", QLL, QLR);
//        readMatrix("../syntheticData/QLMatrix.dat", QLLen, QLCost);
//        RQL(QLCost, QLSeq);
//    }

    for (int i = 0; i < 10; ++i) {
        double btnk;
        VVI seq, L, R;
        VI Lmate, Rmate;
        PII pair;
        VVD cost;
        printf("==================================================================\n");

//        generateSequence("../syntheticData/sequence.dat", len, lambda);
        readSequence("../syntheticData/sequence.dat", seq);
        splitSequence(seq, L, R);
//        generateMatrix("../syntheticData/matrix.dat", L, R);
        readMatrix("../syntheticData/matrix.dat", len, cost);


        printf("Optimal\n");
        btnk = swapChain(cost, L, R, Lmate, Rmate, pair);
        for (int j = 0; j < L.size(); ++j) {
            printf("%4d", j);
        }
        printf("\n");
        for (int j = 0; j < L.size(); ++j) {
            printf("%4d", Lmate[j]);
        }
        printf("\n\n");
        for (int j = 0; j < R.size(); ++j) {
            printf("%4d", j);
        }
        printf("\n");
        for (int j = 0; j < R.size(); ++j) {
            printf("%4d", Rmate[j]);
        }
        printf("\n%.3f < %d, %d >\n\n", btnk, pair.first, pair.second);
//
        printf("Online Greedy\n");
        btnk = onlineGreedy(cost, L, R, Lmate, Rmate, pair);
        for (int j = 0; j < L.size(); ++j) {
            printf("%4d", j);
        }
        printf("\n");
        for (int j = 0; j < L.size(); ++j) {
            printf("%4d", Lmate[j]);
        }
        printf("\n\n");
        for (int j = 0; j < R.size(); ++j) {
            printf("%4d", j);
        }
        printf("\n");
        for (int j = 0; j < R.size(); ++j) {
            printf("%4d", Rmate[j]);
        }
        printf("\n%.4f < %d, %d >\n\n", btnk, pair.first, pair.second);

        printf("CandGreedy\n");
        btnk = CandGreedy(cost, L, R, Lmate, Rmate, pair);
        for (int j = 0; j < L.size(); ++j) {
            printf("%4d", j);
        }
        printf("\n");
        for (int j = 0; j < L.size(); ++j) {
            printf("%4d", Lmate[j]);
        }
        printf("\n\n");
        for (int j = 0; j < R.size(); ++j) {
            printf("%4d", j);
        }
        printf("\n");
        for (int j = 0; j < R.size(); ++j) {
            printf("%4d", Rmate[j]);
        }
        printf("\n%.4f < %d, %d >\n\n", btnk, pair.first, pair.second);
//
        printf("RQL\n");
        btnk = QL(cost, seq, Lmate, Rmate, pair);
        for (int j = 0; j < L.size(); ++j) {
            printf("%4d", j);
        }
        printf("\n");
        for (int j = 0; j < L.size(); ++j) {
            printf("%4d", Lmate[j]);
        }
        printf("\n\n");
        for (int j = 0; j < R.size(); ++j) {
            printf("%4d", j);
        }
        printf("\n");
        for (int j = 0; j < R.size(); ++j) {
            printf("%4d", Rmate[j]);
        }
        printf("\n%.3f < %d, %d >\n\n", btnk, pair.first, pair.second);
//
        printf("==================================================================\n");
    }


    return 0;
}
