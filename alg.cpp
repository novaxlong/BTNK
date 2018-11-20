//
// Created by 王恺忻 on 2018/11/12.
//

#include <queue>
#include <algorithm>
#include "alg.h"
using namespace std;

double Q[MAXSIZE][MAXSIZE][LEN][LEN];

double onlineGreedy(const VVD& cost, VVI& L, VVI& R, VI& Lmate, VI& Rmate, PII& btnkPair) {
    if (L.size() == 0 || R.size() == 0) return 0.;

    VI Lwait, Rwait;        // restore the index of the waiting l/r in the L/R
    double btnk = 0.;
    Lwait.clear();
    Rwait.clear();
    Lmate = VI(L.size(), -1);
    Rmate = VI(R.size(), -1);

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
    return btnk;
}

double CandGreedy(const VVD& cost, VVI& L, VVI& R, VI& Lmate, VI& Rmate, PII& btnkPair) {
    if (L.size() == 0 || R.size() == 0) return 0.;

    VI Lwait, Rwait;
    double btnk = 0., _btnk;
    Lwait.clear();
    Rwait.clear();
    Lmate = VI(L.size(), -1);
    Rmate = VI(R.size(), -1);

    for (int i = 0, j = 0; i < L.size() || j < R.size(); ) {
        if (j == R.size() || (i != L.size() && L[i][2] < R[j][2])) {
            int time = L[i][2];
            for (int k = 0; k < Lwait.size(); ++k) {
                if (L[Lwait[k]][3] == time) {
                    int rindex = 0;
                    double temp = cost[L[Lwait[k]][1]][R[Rwait[0]][1]];
                    for (int p = 0; p < Rwait.size(); ++p) {
                        if (cost[L[Lwait[k]][1]][R[Rwait[p]][1]] >= 0 && cost[L[Lwait[k]][1]][R[Rwait[p]][1]] < temp) {
                            temp = cost[L[Lwait[k]][1]][R[Rwait[p]][1]];
                            rindex = p;
                        }
                    }
                    btnkPair = (btnk > temp) ? btnkPair : make_pair(L[Lwait[k]][1], R[Rwait[rindex]][1]);
                    btnk = (btnk > temp) ? btnk : temp;
                    Lmate[Lwait[k]] = Rwait[rindex];
                    Rmate[Rwait[rindex]] = Lwait[k];
                    Lwait.erase(Lwait.begin() + k);
                    Rwait.erase(Rwait.begin() + rindex);
                    k--;
                }
            }
            for (int k = 0; k < Rwait.size(); ++k) {
                if (R[Rwait[k]][3] == time) {
                    int lindex = 0;
                    double temp = cost[L[Lwait[0]][1]][R[Rwait[k]][1]];
                    for (int p = 0; p < Lwait.size(); ++p) {
                        if (cost[L[Lwait[p]][1]][R[Rwait[k]][1]] >= 0 && cost[L[Lwait[p]][1]][R[Rwait[k]][1]] < temp) {
                            temp = cost[L[Lwait[p]][1]][R[Rwait[k]][1]];
                            lindex = p;
                        }
                    }
                    btnkPair = (btnk > temp) ? btnkPair : make_pair(L[Lwait[lindex]][1], R[Rwait[k]][1]);
                    btnk = (btnk > temp) ? btnk : temp;
                    Lmate[Lwait[lindex]] = Rwait[k];
                    Rmate[Rwait[k]] = Lwait[lindex];
                    Lwait.erase(Lwait.begin() + lindex);
                    Rwait.erase(Rwait.begin() + k);
                    k--;
                }
            }
            if (Rwait.size() == SIZE) {
                int rindex = 0;
                double temp = cost[L[i][1]][R[Rwait[0]][1]];
                for (int k = 0; k < Rwait.size(); ++k) {
                    if (cost[L[i][1]][R[Rwait[k]][1]] >= 0 && cost[L[i][1]][R[Rwait[k]][1]] < temp) {
                        temp = cost[L[i][1]][R[Rwait[k]][1]];
                        rindex = k;
                    }
                }
                btnkPair = (btnk > temp) ? btnkPair : make_pair(L[i][1], R[Rwait[rindex]][1]);
                btnk = (btnk > temp) ? btnk : temp;
                Lmate[i] = Rwait[rindex];
                Rmate[Rwait[rindex]] = i;
                Rwait.erase(Rwait.begin() + rindex);
            }
            else
                Lwait.push_back(i);
            i++;
        }
        else {
            int time = R[j][2];
            for (int k = 0; k < Lwait.size(); ++k) {
                if (L[Lwait[k]][3] == time) {
                    int rindex = 0;
                    double temp = cost[L[Lwait[k]][1]][R[Rwait[0]][1]];
                    for (int p = 0; p < Rwait.size(); ++p) {
                        if (cost[L[Lwait[k]][1]][R[Rwait[p]][1]] >= 0 && cost[L[Lwait[k]][1]][R[Rwait[p]][1]] < temp) {
                            temp = cost[L[Lwait[k]][1]][R[Rwait[p]][1]];
                            rindex = p;
                        }
                    }
                    btnkPair = (btnk > temp) ? btnkPair : make_pair(L[Lwait[k]][1], R[Rwait[rindex]][1]);
                    btnk = (btnk > temp) ? btnk : temp;
                    Lmate[Lwait[k]] = Rwait[rindex];
                    Rmate[Rwait[rindex]] = Lwait[k];
                    Lwait.erase(Lwait.begin() + k);
                    Rwait.erase(Rwait.begin() + rindex);
                    k--;
                }
            }
            for (int k = 0; k < Rwait.size(); ++k) {
                if (R[Rwait[k]][3] == time) {
                    int lindex = 0;
                    double temp = cost[L[Lwait[0]][1]][R[Rwait[k]][1]];
                    for (int p = 0; p < Lwait.size(); ++p) {
                        if (cost[L[Lwait[p]][1]][R[Rwait[k]][1]] >= 0 && cost[L[Lwait[p]][1]][R[Rwait[k]][1]] < temp) {
                            temp = cost[L[Lwait[p]][1]][R[Rwait[k]][1]];
                            lindex = p;
                        }
                    }
                    btnkPair = (btnk > temp) ? btnkPair : make_pair(L[Lwait[lindex]][1], R[Rwait[k]][1]);
                    btnk = (btnk > temp) ? btnk : temp;
                    Lmate[Lwait[lindex]] = Rwait[k];
                    Rmate[Rwait[k]] = Lwait[lindex];
                    Lwait.erase(Lwait.begin() + lindex);
                    Rwait.erase(Rwait.begin() + k);
                    k--;
                }
            }
            if (Lwait.size() == SIZE) {
                int lindex = 0;
                double temp = cost[L[Lwait[0]][1]][R[j][1]];
                for (int k = 0; k < Lwait.size(); ++k) {
                    if (cost[L[Lwait[k]][1]][R[j][1]] >= 0 && cost[L[Lwait[k]][1]][R[j][1]] < temp) {
                        temp = cost[L[Lwait[k]][1]][R[j][1]];
                        lindex = k;
                    }
                }
                btnkPair = (btnk > temp) ? btnkPair : make_pair(L[Lwait[lindex]][1], R[j][1]);
                btnk = (btnk > temp) ? btnk : temp;
                Lmate[Lwait[lindex]] = j;
                Rmate[j] = Lwait[lindex];
                Lwait.erase(Lwait.begin() + lindex);
            }
            else
                Rwait.push_back(j);
            j++;
        }
    }
    VVI LL, RR;
    VI Lm, Rm;
    PII pair;
    for (int i = 0; i < Lwait.size(); ++i) {
        LL.push_back(L[Lwait[i]]);
    }
    for (int i = 0; i < Rwait.size(); ++i) {
        RR.push_back(R[Rwait[i]]);
    }
    _btnk = swapChain(cost, LL, RR, Lm, Rm, pair);
    btnkPair = (btnk > _btnk) ? btnkPair : make_pair(LL[pair.first][1], RR[pair.second][1]);
    btnk = (btnk > _btnk) ? btnk : _btnk;
    for (int i = 0; i < Lm.size(); ++i) {
        Lmate[LL[i][1]] = RR[Lm[i]][1];
        Rmate[RR[Lm[i]][1]] = LL[i][1];
    }
    return btnk;
}

bool BFS(const VVD& cost, VVI& L, VVI& R, VI& Lmate, VI& Rmate, double btnk, PII& btnkPair, VI& chain) {
    if (L.size() == 0 || R.size() == 0) return false;
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
    if (L.size() == 0 || R.size() == 0) return 0.;
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
    return btnk;
}


void initQ() {
    for (int i = 0; i < MAXSIZE; ++i) {
        for (int j = 0; j < MAXSIZE; ++j) {
            for (int k = 0; k < LEN; ++k) {
                for (int p = k; p < LEN; ++p) {
                    Q[i][j][k][p] = INIT * ((double) rand() / (double) RAND_MAX);
                }
            }
        }
    }
}

int tick(const VVI& seq, VVI& L, VVI& R, int start, int t) {
    int time, i;
    for (i = start; i < seq.size() && (time = seq[i][2]) < seq[start][2] + t; ++i) {
        if (seq[i][0] == 0) L.push_back(seq[i]);
        else R.push_back(seq[i]);
        if (L.size() >= MAXSIZE) L.erase(L.begin());
        if (R.size() >= MAXSIZE) R.erase(R.begin());
        for (int j = 0; j < L.size(); ++j) {
            if (time > L[j][3] || L[j][0] == -1) {
                L.erase(L.begin() + j);
                j--;
            }
        }
        for (int j = 0; j < R.size(); ++j) {
            if (time > R[j][3] || R[j][0] == -1) {
                R.erase(R.begin() + j);
                j--;
            }
        }
    }
    return i;
}

void RQL(const VVD& cost, VVI& seq) {
    int l = (int) seq.size();
    double btnk;
    VVI L, R;

    L.clear();
    R.clear();

    int count = 0, i = 0, size_l, size_r, last_l, last_r;
    i = tick(seq, L, R, 0, LOWER_BOUND);
    count = LOWER_BOUND;

    while (i < l) {
        int time = seq[i][2], lt = count;
        size_l = (int) L.size();
        size_r = (int) R.size();
        double temp = Q[size_l][size_r][count - LOWER_BOUND][count - LOWER_BOUND];
        for (int j = count; j <= UPPER_BOUND; ++j) {
            if (Q[size_l][size_r][count - LOWER_BOUND][j - LOWER_BOUND] > temp) {
                temp = Q[size_l][size_r][count - LOWER_BOUND][j - LOWER_BOUND];
                lt = j;
            }
        }
        if (lt == count) {
            VI Lm, Rm;
            PII pair;
            btnk = swapChain(cost, L, R, Lm, Rm, pair);
            for (int j = 0; j < Lm.size(); ++j) {
                if (Lm[j] == -1) continue;
                if (cost[L[j][1]][R[Lm[j]][1]] != -1.) R[Lm[j]][0] = -1;
            }
            for (int j = 0; j < Rm.size(); ++j) {
                if (Rm[j] == -1) continue;
                if (cost[L[Rm[j]][1]][R[j][1]] != -1.) L[Rm[j]][0] = -1;
            }
            for (int j = 0; j < L.size(); ++j) {
                if (time > L[j][3] || L[j][0] == -1) {
                    L.erase(L.begin() + j);
                    j--;
                }
            }
            for (int j = 0; j < R.size(); ++j) {
                if (time > R[j][3] || R[j][0] == -1) {
                    R.erase(R.begin() + j);
                    j--;
                }
            }
            last_l = (int) L.size();
            last_r = (int) R.size();
            i = tick(seq, L, R, i, LOWER_BOUND);
            size_l = (int) L.size();
            size_r = (int) R.size();
            double maxQ = 0.;
            for (int j = 0; j < LEN; ++j) {
                if (Q[size_l][size_r][0][j] > maxQ) {
                    maxQ = Q[size_l][size_r][0][j];
                }
            }
            Q[last_l][last_r][count - LOWER_BOUND][lt - LOWER_BOUND] += ALPHA * (1./btnk + maxQ - Q[last_l][last_r][count - LOWER_BOUND][lt - LOWER_BOUND]);
            count = LOWER_BOUND;
        }
        else {
            last_l = size_l;
            last_r = size_r;
            i = tick(seq, L, R, i, 1);
            size_l = (int) L.size();
            size_r = (int) R.size();
            double maxQ = 0.;
            for (int j = 0; j < LEN; ++j) {
                if (Q[size_l][size_r][count - LOWER_BOUND + 1][j] > maxQ) {
                    maxQ = Q[size_l][size_r][count - LOWER_BOUND + 1][j];
                }
            }
            Q[last_l][last_r][count - LOWER_BOUND][lt - LOWER_BOUND] += ALPHA * (maxQ - Q[last_l][last_r][count - LOWER_BOUND][lt - LOWER_BOUND]);
            count++;
        }
    }
}

double QL(const VVD& cost, VVI& seq, VI& Lmate, VI& Rmate, PII& btnkPair) {
    int n = (int) cost.size();
    int l = (int) seq.size();
    double btnk = 0., _btnk;
    VVI L, R;

    Lmate = VI(n, -1);
    Rmate = VI(n, -1);
    L.clear();
    R.clear();

    int count = 0, i = 0, size_l, size_r, last_l, last_r;
    i = tick(seq, L, R, 0, LOWER_BOUND);
    count = LOWER_BOUND;

    while (i < l) {
        int time = seq[i][2], lt = count;
        size_l = (int) L.size();
        size_r = (int) R.size();
        double temp = Q[size_l][size_r][count - LOWER_BOUND][count - LOWER_BOUND];
        for (int j = count; j <= UPPER_BOUND; ++j) {
            if (Q[size_l][size_r][count - LOWER_BOUND][j - LOWER_BOUND] > temp) {
                temp = Q[size_l][size_r][count - LOWER_BOUND][j - LOWER_BOUND];
                lt = j;
            }
        }
        if (lt == count) {
            VI Lm, Rm;
            PII pair;
            _btnk = swapChain(cost, L, R, Lm, Rm, pair);
            btnkPair = (btnk > _btnk) ? btnkPair : make_pair(L[pair.first][1], R[pair.second][1]);
            btnk = (btnk > _btnk) ? btnk : _btnk;
            for (int j = 0; j < Lm.size(); ++j) {
                if (Lm[j] == -1) continue;
                if (cost[L[j][1]][R[Lm[j]][1]] != -1.)  {
                    Lmate[L[j][1]] = R[Lm[j]][1];
                    Rmate[R[Lm[j]][1]] = L[j][1];
                    R[Lm[j]][0] = -1;
                }
            }
            for (int j = 0; j < Rm.size(); ++j) {
                if (Rm[j] == -1) continue;
                if (cost[L[Rm[j]][1]][R[j][1]] != -1.) {
                    Lmate[L[Rm[j]][1]] = R[j][1];
                    Rmate[R[j][1]] = L[Rm[j]][1];
                    L[Rm[j]][0] = -1;
                }
            }
            for (int j = 0; j < L.size(); ++j) {
                if (time > L[j][3] || L[j][0] == -1) {
                    L.erase(L.begin() + j);
                    j--;
                }
            }
            for (int j = 0; j < R.size(); ++j) {
                if (time > R[j][3] || R[j][0] == -1) {
                    R.erase(R.begin() + j);
                    j--;
                }
            }
            i = tick(seq, L, R, i, LOWER_BOUND);
            count = LOWER_BOUND;
        }
        else {
            i = tick(seq, L, R, i, 1);
            count++;
        }
    }
    VI Lm, Rm;
    PII pair;
    _btnk = swapChain(cost, L, R, Lm, Rm, pair);
    btnkPair = (btnk > _btnk) ? btnkPair : make_pair(L[pair.first][1], R[pair.second][1]);
    btnk = (btnk > _btnk) ? btnk : _btnk;
    for (int j = 0; j < Lm.size(); ++j) {
        if (Lm[j] == -1) continue;
        if (cost[L[j][1]][R[Lm[j]][1]] != -1.)  {
            Lmate[L[j][1]] = R[Lm[j]][1];
            Rmate[R[Lm[j]][1]] = L[j][1];
            R[Lm[j]][0] = -1;
        }
    }
    for (int j = 0; j < Rm.size(); ++j) {
        if (Rm[j] == -1) continue;
        if (cost[L[Rm[j]][1]][R[j][1]] != -1.) {
            Lmate[L[Rm[j]][1]] = R[j][1];
            Rmate[R[j][1]] = L[Rm[j]][1];
            L[Rm[j]][0] = -1;
        }
    }
    return btnk;
}
