// poj 1062
#include <iostream>
#include <vector>
#include <climits>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <set>
#define S(fmt, ...) do {if(true) printf(fmt"\n", ##__VA_ARGS__);} while(false)
using namespace std;
const int MAXN = 110;
const int INF = 1<<28;
#define ABS(x) ((x)>0?(x):-(x))
#define MIN(a, b) (a) < (b) ? (a) : (b)
#define MAX(a, b) (a) > (b) ? (a) : (b)
int M, N;
struct node {
    int price, level, x;
    vector< pair<int, int> > tx;
};
node A[MAXN];

int matrix[MAXN][MAXN];
int mark[MAXN], dist[MAXN];

void buildG() {
    int i, j;
    for (i = 0; i < MAXN; i++)
    for (j = 0; j < MAXN; j++)
        matrix[i][j] = INF;
    for (i = 1; i <= N; ++i) {
        node& no = A[i];
        matrix[N + 1][i] = no.price;
        for (j = 0; j < no.x; j++) {
            int idx = no.tx[j].first;
            node& no2 = A[idx];
            int p = no.tx[j].second;
            if (ABS(no.level - no2.level) <= M) {
                matrix[idx][i] = p;
            }
        }
    }
}
int Dijkstra(int f) {
    int i, k, min;
    int x = N + 1; int y = 1;
    memset(mark, 0, sizeof(mark));
    for (i = 1; i <= N + 1; i++) dist[i] = matrix[x][i];
    mark[x] = 1;
    do {
        min = INF;
        k = 0;
        for (i = 1; i <= N + 1; i++) if (mark[i] == 0 && dist[i] < min) {
            min = dist[i];
            k = i;
        }
        if (k) {
            mark[k] = 1;
            for (i = 1; i <= N + 1;i++) if (matrix[k][i] != INF &&
                                            (min + matrix[k][i] < dist[i]) &&
                                            A[i].level - f <= M && A[i].level - f >= 0 &&
                                            A[k].level - f <= M && A[k].level - f >= 0
                                            ) {

                dist[i] = min + matrix[k][i];
            }
        }
    } while(k);
    return dist[y];
}

int main()
{
    ios::sync_with_stdio(false);
    while(cin>>M>>N) {
        set<int> sl; sl.clear();
        vector<int> vr; vr.clear();

        for (int i = 1; i <= N; ++i) {
            cin>>A[i].price>>A[i].level>>A[i].x;
            if (sl.find(A[i].level) == sl.end()) sl.insert(A[i].level);
            for (int j = 0; j < A[i].x; ++j) {
                int t, v;
                cin>>t>>v;
                A[i].tx.push_back(make_pair(t, v));
            }
        }
        set<int>::iterator it = sl.begin();
        buildG();
        while(it != sl.end()) {
            vr.push_back(Dijkstra(*it));
            it++;
        }
        cout<<*min_element(vr.begin(), vr.end())<<endl;
    }
    return 0;
}

