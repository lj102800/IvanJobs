#include <cstdio>
#include <cstdlib>

int main() {
	int n, m;
	while (scanf("%d%d", &n, &m) != EOF) {
		int i, j;
		int val, ri, rj, rv;
		ri = rj = -1;
		rv = 0;
		for (i = 0; i < n; ++i) {
			for (j = 0; j < m; ++j) {
				scanf("%d", &val);
				if (abs(val) > abs(rv)) {
					rv = val;
					ri = i;
					rj = j;
				}	
			}
		}
		printf("%d %d %d\n", ri + 1, rj + 1, rv);
	}
	return 0;
}
