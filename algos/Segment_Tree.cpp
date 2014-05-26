/**
 * Segment Tree code template from boj.
 * Tree build on an MAX segment(or range). The leaf node can be [i, i] or [i , i + 1] depends on u.
 * and every node is an range, and it's two sons and split by, [a, (a + b)/2] [(a + b)/2 + 1, b], if the 
 * current node is [a, b];
 * or split by [a, (a + b)/2] , [(a + b)/2 , b], when stop on range length == 1 
 */
struct node {
	int l, r;
	int c;
}T[3*MAXN];

void build(int l, int r, int k) {
	if (l == r) {
		T[k].l = l; T[k].r = r; T[k].c = 0; return ;
	}
	int m = (l + r)>>1;
	T[k].l = l; T[k].r = r; T[k].c = 0;
	build(l, m, 2*k);
	build(m + 1, r, 2*k + 1);
}

void insert(int d, int k) {
	if (T[k].l == T[k].r && d == T[k].l) {
		T[k].c += 1;
		return ;
	}
	int m = (T[k].l + T[k].r)>>1;
	if (d <= m) insert(d, 2 * k);
	else insert(d, 2 * k + 1);
	T[k].c = T[2 * k].c + T[2 * k + 1].c;
}

void search(int d, int k, int& ans) {
	if (T[k].l == T[k].r) {
		ans = T[k].l;
		return ;
	}
	int m = (T[k].l + T[k].r) >> 1;
	if (d > T[2 * k].c)
		search(d - T[2 * k].c, 2 * k + 1, ans);
	else search(d, 2 * k, ans);
}


