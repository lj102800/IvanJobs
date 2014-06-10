// hdu 1102
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>

#define REP(i, n) for(int _n = n, i = 0; i < _n; i++)
#define FOR(i, a, b) for(int i = (a), _b = (b); i <= _b; i++)
#define Max(a, b) ((a) > (b) ? (a) : (b))
#define Min(a, b) ((a) < (b) ? (a) : (b))
#define L(fmt, ...) do {if(false) printf(fmt"\n", ##__VA_ARGS__);} while(false)

struct EDGE {
	int bv, tv;
	int w;
};

EDGE edges[10100];
int matrix[110][110];
bool built[110][110];
int UF[110], V, E;
int sum;
int root(int v) {
    if (UF[v] != v) {
        UF[v] = root(UF[v]);
    }
    return UF[v];
}

int UDlesser(const void* a, const void* b) {
	EDGE* pa = (EDGE*) a;
	EDGE* pb = (EDGE*) b;
	if (pa->w < pb->w) return -1;
	else if (pa->w > pb->w) return 1;
	else return 0;
}

void kruskal(int cnt){
	int v1, v2, i, j;
	sum = i = j = 0;
	L("cnt:%d, E:%d", cnt, E);
	while(j < cnt && i < E) {
		v1 = root(edges[i].bv);
		v2 = root(edges[i].tv);
		L("v1:%d, v2:%d", v1, v2);
		if (v1 != v2) {
			sum += edges[i].w;
			UF[v1] = v2;
			j++;
		}
		i++;
	}
}

int main() {
	while(scanf("%d", &V) != EOF) {
	    FOR(i, 1, V) UF[i] = i;
        FOR(i, 1, V) {
            FOR(j, 1, V) {
                scanf("%d", &matrix[i][j]);
            }
        }
        memset(built, false, sizeof(built));
        int Q; int cnt = 0;
        scanf("%d", &Q);
        REP(i, Q) {
            int a, b;
            scanf("%d%d", &a, &b);
            built[a][b] = built[b][a] = true;
            int ra = root(a);
            int rb = root(b);
            if (ra != rb) {
                cnt++;
                UF[ra] = rb;
            }
        }
        E = 0;
        FOR(i, 1, V) {
            FOR(j, i + 1, V) {
                if (!built[i][j]) {
                    edges[E].bv = i;
                    edges[E].tv = j;
                    edges[E].w = matrix[i][j];
                    E++;
                }
            }
        }
        qsort(edges, E, sizeof(EDGE), UDlesser);
        kruskal(V - cnt);
        printf("%d\n", sum);

	}
	return 0;
}
