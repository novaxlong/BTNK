//
// Created by 王恺忻 on 2018/11/10.
//

#include <random>
#include "data.h"

using namespace std;

void generateSequence(const char *fileName, int len, int lambda) {
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
    int taskId = 0, workerId = 0, start = 0, end = 0, duration = 0;
    random_device rd;
    mt19937 gen(rd());
    poisson_distribution<> pd(lambda);
    uniform_int_distribution<> uid(LOWER_BOUND_SEQ, UPPER_BOUND_SEQ);

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
            duration = uid(gen);
            end = start + duration;
            fprintf(fp, "0 %d %d %d\n", taskId++, start, end);
        }
        for (int j = 0; j < worker_k; ++j) {
            duration = uid(gen);
            end = start + duration;
            fprintf(fp, "1 %d %d %d\n", workerId++, start, end);
        }
        taskLen += task_k;
        workerLen += worker_k;
        start++;
        num += k;
    }
    fclose(fp);
}

void readSequence(const char *fileName, VVI &seq) {
    FILE *fp = fopen(fileName, "r");
    if (fp == nullptr) {
        printf("Cannot open file %s\n.", fileName);
        exit(-1);
    }
    VI temp;
    int type, id, start, end, matched_time = 0;
    while (fscanf(fp, "%d %d %d %d", &type, &id, &start, &end) == 4) {
        temp.push_back(type);           // 0
        temp.push_back(id);             // 1
        temp.push_back(start);          // 2
        temp.push_back(end);            // 3
        temp.push_back(matched_time);   // 4
        seq.push_back(temp);
        temp.clear();
    }
    fclose(fp);
}

void splitSequence(VVI &seq, VVI &L, VVI &R) {
    L.clear();
    R.clear();
    for (int i = 0; i < seq.size(); ++i) {
        if (seq[i][0] == 0)
            L.push_back(seq[i]);
        else
            R.push_back(seq[i]);
    }
}

/*
 * Modified on 2018.11.22 by Wang.
 * TODO: To generate the distance between node_l and node_r.
 */
void generateMatrix(const char *fileName, VVI &L, VVI &R) {
    FILE *fp = fopen(fileName, "w");
    if (fp == nullptr) {
        printf("Cannot open file %s\n.", fileName);
        exit(-1);
    }
    double cost;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> urd(LOWER_BOUND_DIS, UPPER_BOUND_DIS);

    for (int i = 0; i < L.size(); ++i) {
        for (int j = 0; j < R.size(); ++j) {
            if ((L[i][2] - R[j][3]) * (L[i][3] - R[j][2]) <= 0) cost = urd(gen);
            else cost = -1;
            fprintf(fp, "%7.3lf ", cost);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

/*
 * Modified on 2018.11.22 by Wang.
 * TODO: To read the distance between node_l and node_r.
 */

void readMatrix(const char *fileName, int len, VVD &mat) {
    FILE *fp = fopen(fileName, "r");
    if (fp == nullptr) {
        printf("Cannot open file %s\n.", fileName);
        exit(-1);
    }
    double value;
    mat = VVD((unsigned long) len, VD((unsigned long) len));
    for (int i = 0; i < len; ++i) {
        for (int j = 0; j < len; ++j) {
            fscanf(fp, "%lf ", &value);
            mat[i][j] = value;
        }
    }
    fclose(fp);
}
