// poj 1125
#include <cstdio>
#include <queue>
#include <cstring>
#include <algorithm>
#include <vector>
using namespace std;

#define REP(i, n) for(int _n = n, i = 0; i < _n; i++)
#define FOR(i, a, b) for(int i = (a), _b = (b); i <= _b; i++)
#define Max(a, b) ((a) > (b) ? (a) : (b))
#define Min(a, b) ((a) < (b) ? (a) : (b))
#define L(fmt, ...) do {if(false) printf(fmt"\n", ##__VA_ARGS__);} while(false)
const int INF = 1000000;
int matrix[110][110];
int dist[110][110];
int N;
bool vis[110];
void bfs() {
   queue<int> q;
   memset(vis, false, sizeof(vis));
   q.push(1);
   while(!q.empty()) {
        int t = q.front();
        q.pop();
        vis[t] = true;
        FOR(i, 1, N) {
            if (!vis[i] && (matrix[t][i] < INF || matrix[i][t] < INF)) {
                q.push(i);
            }
        }
   }
}
bool connected() {
    bfs();
    FOR(i, 1, N) {
        if (!vis[i]) return false;
    }
    return true;
}

void floyd() {
    FOR(i, 1, N) {
        FOR(j, 1, N) {
            if(i != j)
                dist[i][j] = matrix[i][j];
            else dist[i][j] = 0;
        }
    }
    FOR(k, 1, N) {
        FOR(i, 1, N) {
            FOR(j, 1, N) {
                if (dist[i][j] > dist[i][k] + dist[k][j]) {
                    dist[i][j] = dist[i][k] + dist[k][j];
                }
            }
        }
    }
}

int main() {
    while(scanf("%d", &N) && N) {
        REP(i, 110) {
            REP(j, 110) {
                matrix[i][j] = INF;
            }
        }
        int n;
        FOR(i, 1, N) {
            scanf("%d", &n);
            REP(j, n) {
                int b, t;
                scanf("%d%d", &b, &t);
                matrix[i][b] = t;
            }
        }
        if (!connected()) {
            printf("disjoint\n");
        } else {
            floyd();
            vector<int> vr; vr.clear();
            FOR(i, 1, N) {
                vr.push_back(*max_element(dist[i] + 1, dist[i] + N + 1));
            }
            int midx = 0; int mv = INF;
            REP(i, vr.size()) {
                if (vr[i] < mv) {
                    mv = vr[i];
                    midx = i;
                }
            }
            printf("%d %d\n", midx + 1, vr[midx]);
        }
    }
    return 0;
}
