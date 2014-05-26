/**
 * Segment Tree code template from boj.
 * Tree build on an MAX segment(or range). The leaf node can be [i, i] or [i , i + 1] depends on u.
 * and every node is an range, and it's two sons and split by, [a, (a + b)/2] [(a + b)/2 + 1, b], if the
 * current node is [a, b];
 * or split by [a, (a + b)/2] , [(a + b)/2 , b], when stop on range length == 1
 */
#include <cstdio>
#define MAXN 200100
int MAX(int a, int b) {return (a) > (b) ? (a) : (b);}
#define L(fmt, ...) do {if(false) printf(fmt"\n", ##__VA_ARGS__);} while(false)
struct node {
	int l, r;
	int c;
}T[3*MAXN];
int A[MAXN];
void build(int l, int r, int k) {
    L("build(%d, %d, %d)", l, r, k);
	if (l == r) {
		T[k].l = l; T[k].r = r; T[k].c = A[T[k].l]; return ;
	}
	int m = (l + r)>>1;
	T[k].l = l; T[k].r = r; T[k].c = 0;
	build(l, m, 2*k);
	build(m + 1, r, 2*k + 1);
	T[k].l = l; T[k].r = r;
	T[k].c = MAX(T[2 * k].c , T[2 * k + 1].c);
}

void update(int d, int v, int k) {
    L("update(%d, %d, %d)", d, v, k);
	if (T[k].l == T[k].r ) {
	    if (d == T[k].l)
            T[k].c = T[k].c < v ? v : T[k].c;
		return ;
	}
	int m = (T[k].l + T[k].r)>>1;
	if (d <= m) {
	    update(d, v, 2 * k);
	    T[k].c = T[2 * k].c > T[k].c ? T[2 * k].c : T[k].c;
	}
	else if(d >= m + 1){
	    update(d, v, 2 * k + 1);
	    T[k].c = T[2 * k + 1].c > T[k].c ? T[2 * k + 1].c : T[k].c;
	}
}

int query(int l, int r,  int k) {
    L("query(%d, %d, %d)", l, r, k);
	if (T[k].l == l &&  T[k].r == r) {
		return T[k].c;
	}
	int m = (T[k].l + T[k].r) >> 1;
	if (m >= r) {
	    return query(l, r, 2 * k);
	} else if (m + 1 <= l) {
        return query(l, r, 2 * k + 1);
	} else {
        return MAX(query(l, m, 2 * k), query(m + 1, r, 2 * k + 1));
	}
}

int main() {
    int N, M;
    while(scanf("%d%d", &N, &M) != EOF) {
        int i;
        for (i = 1; i <= N; ++i) scanf("%d", &A[i]);
        build(1, N, 1);
        for (i = 0; i < M; ++i) {
            char cmd[2];
            int A, B;
            scanf("%s", cmd);
            scanf("%d%d", &A, &B);
            switch(cmd[0]) {
                case 'Q':
                    printf("%d\n", query(A, B, 1));
                    break;
                case 'U':
                    update(A, B, 1);
                    break;
            }
        }
    }
    return 0;
}

