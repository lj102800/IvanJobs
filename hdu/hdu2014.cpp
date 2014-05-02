#include <cstdio>
#include <algorithm>

using namespace std;
int mat[150];

int main() {
	int n;
	while(scanf("%d", &n) != EOF) {
		int i;
		for (i = 0; i < n; ++i) {
			scanf("%d", &(mat[i]));
		}
		sort(mat, mat + n);
		int s = 0;
		for (i = 1; i <= n - 2; ++i) {
			s += mat[i];
		}
		printf("%.2f\n", (float)s / (n - 2));
	}
	return 0;
}
