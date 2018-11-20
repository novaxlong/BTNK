//
// Created by 王恺忻 on 2018/11/12.
//

#include <queue>
#include <algorithm>
#include "alg.h"

using namespace std;

bool cmp(const VD& i, const VD& j) {
    return (i[0] < j[0]);
}

double onlineGreedy(const VVD& cost, VVI& L, VVI& R, VI& Lmate, VI& Rmate, PII& btnkPair) {
    VI Lwait, Rwait;        // restore the index of the waiting l/r in the L/R
    double btnk = 0.;
    Lwait.clear();
    Rwait.clear();
    Lmate = VI(L.size(), -1);
    Rmate = VI(R.size(), -1);
    btnkPair = make_pair(-1, -1);
    for (int i = 0, j = 0; i < L.size() || j < R.size(); ) {
        if (j == R.size() || (i != L.size() && L[i][2] < R[j][2])) {
            for (int k = 0; k < Lwait.size(); ++k) {
                if (L[Lwait[k]][3] < L[i][2]) {
                    printf("Discard task %d.\n", L[Lwait[k]][1]);
                    Lwait.erase(Lwait.begin() + k);
                    k--;
                }
            }
            for (int k = 0; k < Rwait.size(); ++k) {
                if (R[Rwait[k]][3] < L[i][2]) {
                    printf("Discard worker %d.\n", R[Rwait[k]][1]);
                    Rwait.erase(Rwait.begin() + k);
                    k--;
                }
            }
            if (!Rwait.empty()) {
                double minCost = 1.0;
                int minIndex = -1, rIndex, kIndex = -1;
                for (int k = 0; k < Rwait.size(); ++k) {
                    rIndex = Rwait[k];
                    if (cost[L[i][1]][R[rIndex][1]] >= 0 && cost[L[i][1]][R[rIndex][1]] <= minCost) {
                        minCost = cost[L[i][1]][R[rIndex][1]];
                        minIndex = rIndex;
                        kIndex = k;
                    }
                }
                Lmate[i] = minIndex;
                Rmate[minIndex] = i;
                btnkPair = (btnk > minCost) ? btnkPair : make_pair(i, minIndex);
                btnk = (btnk > minCost) ? btnk : minCost;
                Rwait.erase(Rwait.begin() + kIndex);
            }
            else
                Lwait.push_back(i);
            i++;
        }
        else {
            for (int k = 0; k < Lwait.size(); ++k) {
                if (L[Lwait[k]][3] < R[j][2]) {
                    printf("Discard task %d.\n", L[Lwait[k]][1]);
                    Lwait.erase(Lwait.begin() + k);
                    k--;
                }
            }
            for (int k = 0; k < Rwait.size(); ++k) {
                if (R[Rwait[k]][3] < R[j][2]) {
                    printf("Discard worker %d.\n", R[Rwait[k]][1]);
                    Rwait.erase(Rwait.begin() + k);
                    k--;
                }
            }
            if (!Lwait.empty()) {
                double minCost = 1.0;
                int minIndex = -1, lIndex, kIndex = -1;
                for (int k = 0; k < Lwait.size(); ++k) {
                    lIndex = Lwait[k];
                    if (cost[L[lIndex][1]][R[j][1]] >= 0 && cost[L[lIndex][1]][R[j][1]] <= minCost) {
                        minCost = cost[L[lIndex][1]][R[j][1]];
                        minIndex = lIndex;
                        kIndex = k;
                    }
                }
                Lmate[minIndex] = j;
                Rmate[j] = minIndex;
                btnkPair = (btnk > minCost) ? btnkPair : make_pair(minIndex, j);
                btnk = (btnk > minCost) ? btnk : minCost;
                Lwait.erase(Lwait.begin() + kIndex);
            }
            else
                Rwait.push_back(j);
            j++;
        }
    }
    for (int i = 0; i < L.size(); ++i) {
        printf("%4d", L[i][1]);
    }
    printf("\n");
    for (int i = 0; i < Lmate.size(); ++i) {
        if (Lmate[i] == -1)
            printf("%4d", -1);
        else
            printf("%4d", R[Lmate[i]][1]);
    }
    printf("\n\n");
    for (int i = 0; i < R.size(); ++i) {
        printf("%4d", R[i][1]);
    }
    printf("\n");
    for (int i = 0; i < Rmate.size(); ++i) {
        if (Rmate[i] == -1)
            printf("%4d", -1);
        else
            printf("%4d", L[Rmate[i]][1]);
    }
    printf("\n\n");
    return btnk;
}

bool BFS(const VVD& cost, VVI& L, VVI& R, VI& Lmate, VI& Rmate, double btnk, PII& btnkPair, VI& chain) {
    int l = btnkPair.first, r = btnkPair.second;
    VVI queue;
    bool update = false;
    VI dad(R.size());
    vector<bool> seen = vector<bool>(R.size(), false);

    VI temp = {0, btnkPair.first};
    queue.push_back(temp);
    Lmate[l] = -1;
    Rmate[r] = -1;

    while (!queue.empty()) {
        VI head = queue[0];
        if (head[0] == 0) {
            for (int i = 0; i < R.size(); ++i) {
                if (!seen[i] && cost[L[head[1]][1]][R[i][1]] >= 0 && cost[L[head[1]][1]][R[i][1]] < btnk) {
                    dad[i] = head[1];
                    seen[i] = true;
                    temp = {1, i};
                    queue.push_back(temp);
                }
            }
        }
        else {
            if (Rmate[head[1]] == -1) {
                r = head[1];
                update = true;
                break;
            }
            temp = {0, Rmate[head[1]]};
            queue.push_back(temp);
        }
        queue.erase(queue.begin());
    }
    if (update) {
        l = dad[r];
        chain.push_back(r);
        chain.push_back(l);
        while (l != btnkPair.first) {
            r = Lmate[l];
            l = dad[r];
            chain.push_back(r);
            chain.push_back(l);
        }
        return true;
    }
    Lmate[l] = r;
    Rmate[r] = l;
    return false;
}

double swapChain(const VVD& cost, VVI& L, VVI& R, VI& Lmate, VI& Rmate, PII& btnkPair) {
    double btnk = onlineGreedy(cost, L, R, Lmate, Rmate, btnkPair);
    VI chain;
    while (BFS(cost, L, R, Lmate, Rmate, btnk, btnkPair, chain)) {
        for (int i = 0; i < chain.size(); ++i) {
            if (i & 1) {
                Lmate[chain[i]] = chain[i-1];
                Rmate[chain[i-1]] = chain[i];
            }
        }
        btnk = 0.;
        for (int i = 0; i < Lmate.size(); ++i) {
            if (Lmate[i] == -1) continue;
            btnkPair = (btnk > cost[L[i][1]][R[Lmate[i]][1]] ? btnkPair : make_pair(i, Lmate[i]));
            btnk = (btnk > cost[L[i][1]][R[Lmate[i]][1]] ? btnk : cost[L[i][1]][R[Lmate[i]][1]]);
        }
        chain.clear();
    }
    for (int i = 0; i < L.size(); ++i) {
        printf("%4d", L[i][1]);
    }
    printf("\n");
    for (int i = 0; i < Lmate.size(); ++i) {
        if (Lmate[i] == -1)
            printf("%4d", -1);
        else
            printf("%4d", R[Lmate[i]][1]);
    }
    printf("\n\n");
    for (int i = 0; i < R.size(); ++i) {
        printf("%4d", R[i][1]);
    }
    printf("\n");
    for (int i = 0; i < Rmate.size(); ++i) {
        if (Rmate[i] == -1)
            printf("%4d", -1);
        else
            printf("%4d", L[Rmate[i]][1]);
    }
    printf("\n");
    return btnk;
}




