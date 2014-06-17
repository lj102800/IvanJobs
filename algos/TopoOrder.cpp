// hdu 2647
#include <cstdio>
#include <cstring>
#include <vector>
#include <queue>
using namespace std;

#define REP(i, n) for(int _n = n, i = 0; i < _n; i++)
#define FOR(i, a, b) for(int i = (a), _b = (b); i <= _b; i++)
#define Max(a, b) ((a) > (b) ? (a) : (b))
#define Min(a, b) ((a) < (b) ? (a) : (b))
#define L(fmt, ...) do {if(false) printf(fmt"\n", ##__VA_ARGS__);} while(false)
#define MAXN 10010
int n, m;
int cnt[MAXN], value[MAXN];
vector<int> adj[MAXN];
int TopoOrder() {
    queue<int> q;
    int i;
    int res = 0; int rn = 0;
    for (i = 0; i < n; ++i) if (cnt[i] == 0) {
        L("i:%d", i);
        value[i] = 888;
        q.push(i);
    }
    while(!q.empty()) {
        int j = q.front();
        q.pop(); rn++;
        vector<int>& si = adj[j];
        vector<int>::iterator it = si.begin();
        while(it != si.end()) {
            int k = *it;
            if ((--cnt[k]) == 0) {
                value[k] = value[j] + 1;
                q.push(k);
            }
            it++;
        }
    }
    if (rn != n) return -1;
    REP(i, n) {
        res += value[i];
    }
    return res;
}

int main() {
    while(scanf("%d%d", &n, &m) != EOF ) {
        memset(cnt, 0, sizeof(cnt));
        memset(value, 0, sizeof(value));
        REP(i, MAXN) adj[i].clear();
        REP(i, m) {
            int a, b;
            scanf("%d%d", &a, &b);
            a--; b--;
            adj[b].push_back(a);
            cnt[a]++;
        }
        int ret = TopoOrder();
        printf("%d\n", ret);
    }
    return 0;
}
