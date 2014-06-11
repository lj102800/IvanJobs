// hdu 4081
#include <cstdio>
#include <climits>
#include <cstring>
#include <cmath>
#include <vector>
using namespace std;

#define REP(i, n) for(int _n = n, i = 0; i < _n; i++)
#define FOR(i, a, b) for(int i = (a), _b = (b); i <= _b; i++)
#define Max(a, b) ((a) > (b) ? (a) : (b))
#define Min(a, b) ((a) < (b) ? (a) : (b))
#define SSquare(x1, y1, x2, y2) (((x1)-(x2))*((x1)-(x2))+((y1)-(y2))*((y1)-(y2)))
#define L(fmt, ...) do {if(false) printf(fmt"\n", ##__VA_ARGS__);} while(false)
#define MAXV 1010

const double eps = 0.0000001;
int cost[MAXV][MAXV]; // distance between each city pair, in sum square
int X[MAXV], Y[MAXV], P[MAXV]; // population for each city
bool vis[MAXV];
int Father[MAXV];
bool MST[MAXV][MAXV];
int Mx[MAXV][MAXV];
int lowc[MAXV];
vector<double> out;
int n; // number of cities
double prim() {
    int p;
    int minc;
    double res = 0.0;
    memset(vis, false, sizeof(vis));
    memset(MST, false, sizeof(MST));
    memset(Mx, 0, sizeof(Mx));
    vis[0] = true;
    FOR(i, 1, n - 1) {
        lowc[i] = cost[0][i];
        Father[i] = 0;
    }
    FOR(i, 1, n - 1) {
        minc = INT_MAX; p = -1;
        REP(j, n) if (!vis[j] && lowc[j] < minc) {
            minc = lowc[j];
            p = j;
        }
        L("find minc:%d, p:%d, ", minc, p);
        res += sqrt((double)minc);
        MST[p][Father[p]] = MST[Father[p]][p] = true;
        vis[p] = true;
        REP(j, n) {
            if (vis[j] && j != p) {
                Mx[p][j] = Mx[j][p] = Max(Mx[j][Father[p]], lowc[p]);
            }
            if (!vis[j] && lowc[j] > cost[p][j]) {
                lowc[j] = cost[p][j];
                Father[j] = p;
            }
        }
    }
    return res;
}

void makeMx() {
    memset(Mx, 0, sizeof(Mx));
    REP(i, n) {
        int j = i; int m = -1.0;
        while(Father[j] != j) {
            if (cost[Father[j]][j] > m) {
                m = cost[Father[j]][j];
            }
            Mx[i][Father[j]] = Mx[Father[j]][i] = m;
            j = Father[j];
        }
    }
}

int main() {
    int t;
    scanf("%d", &t);
    while(t--) {
        scanf("%d", &n);
        REP(i, n) {
            scanf("%d%d%d", &X[i], &Y[i], &P[i]);
        }
        REP(i, n) {
            FOR(j, i + 1, n - 1) {
                cost[i][j] = cost[j][i] = SSquare(X[i], Y[i], X[j], Y[j]);
                L("cost[%d][%d]:%d", i, j, cost[i][j]);
            }
        }
        double res = prim();
        L("res:%lf", res);
        out.clear();
        REP(i, n) {
            FOR(j, i + 1, n - 1) {
                int pa = P[i] + P[j];
                L("pa:%d", pa);
                if (MST[i][j]) {
                    L("%d -> %d in MST", i, j);
                    L("cost[%d][%d]:%d", i, j, cost[i][j]);
                    double r1 = (double)pa / (res - sqrt((double)cost[i][j]));
                    L("r1:%lf", r1);
                    out.push_back(r1);
                } else {
                    L("%d -> %d not in MST", i, j);
                    L("Mx[%d][%d]:%d", i, j, Mx[i][j]);
                    double r2 = (double)pa / (res - sqrt((double)Mx[i][j]));
                    L("r2:%lf", r2);
                    out.push_back(r2);
                }
            }
        }
        int len = out.size();
        double r = -1.0;
        REP(i, len) if (out[i] - r > eps) {
            r = out[i];
        }
        printf("%.2lf\n", r);
    }
    return 0;
}

