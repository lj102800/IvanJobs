#include <cstdio>

int mat[110000];

//segment seg[110000];

struct segment {
	int l, r;
};

struct segment seg[110000];
int main() {
	int T;
	int n;
	scanf("%d", &T);
	int cid;
	for (cid = 1; cid <= T; ++cid) {
		scanf("%d", &n);
		int i;
		for (i = 1; i <= n; ++i) {
			scanf("%d", &(mat[i]));
		}
		
		seg[1].l = 1;
		seg[1].r = 1;
		for (i = 2; i <= n; ++i) {
			if (mat[i - 1] + mat[i] > mat[i]) {
				mat[i] = mat[i - 1] + mat[i];
				seg[i].l = seg[i - 1].l;
				seg[i].r = i;	
			} else {
				mat[i] = mat[i];
				seg[i].l = seg[i].r = i;
			}
		}
		int maxR = mat[1];
		int maxId = 1;	
		for (i = 2 ; i <= n; ++i) {
			if (mat[i] > maxR) {
				maxR = mat[i];
				maxId = i;
			}
		}
		if (cid != 1) {
			printf("\n");
		}
		printf("Case %d:\n", cid);
		printf("%d %d %d\n", maxR, seg[maxId].l, seg[maxId].r);		
	}

	return 0;
}
