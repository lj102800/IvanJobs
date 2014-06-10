// poj 1860
#include <cstdio>
#include <cstring>
#include <vector>
using namespace std;
#define REP(i, n) for(int _n = n, i = 0; i < _n; i++)
#define FOR(i, a, b) for(int i = (a), _b = (b); i <= _b; i++)
#define Max(a, b) ((a) > (b) ? (a) : (b))
#define Min(a, b) ((a) < (b) ? (a) : (b))
#define L(fmt, ...) do {if(true) printf(fmt"\n", ##__VA_ARGS__);} while(false)
const double INF = -100.0;
const double eps = 0.0000001;
struct EDGE {
    EDGE():r(INF), c(INF) {}
    EDGE(int _A, int _B, double _r, double _c) : A(_A), B(_B), r(_r), c(_c) {}
    int A, B;
    double r, c;
};
int N, M, S;
double V;

vector<EDGE> edges;
double dist[110];


bool Bellman() {
    FOR(i, 1, N) {
        dist[i] = INF;
    }
    dist[S] = V;
    while(true) {
        int len = edges.size();
        bool bgoing = false;
        REP(i, len) {
            EDGE& ed = edges[i];
            int a, b; double r, c;
            a = ed.A; b = ed.B; r = ed.r; c = ed.c;
            if (dist[a] - c > eps) {
                if (dist[b] < (dist[a] - c) * r) {
                    dist[b] = (dist[a] - c) * r;
                    bgoing = true;
                    if (dist[S] - V > eps) return true;
                }
            }
        }
        if (!bgoing) break;
    }
    return false;
}

int main() {
    while(scanf("%d%d%d%lf", &N, &M, &S, &V) != EOF){
        edges.clear();
        REP(i, M) {
            int a, b; double Rab, Rba, Cab, Cba;
            scanf("%d%d%lf%lf%lf%lf", &a, &b, &Rab, &Cab, &Rba, &Cba);
            edges.push_back(EDGE(a, b, Rab, Cab));
            edges.push_back(EDGE(b, a, Rba, Cba));
        }
        bool bok = Bellman();
        if (bok) printf("YES\n");
        else printf("NO\n");
    }
    return 0;
}
