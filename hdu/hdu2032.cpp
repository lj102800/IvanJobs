#include <cstdio>
int a[40][40];

void make() {
	a[1][1] = a[2][1] = a[2][2] = 1;
	int i, j;
	for (i = 3; i <= 30; ++i) {
		a[i][1] = 1;
		a[i][i] = 1;
		for (j = 2; j < i; ++j) {
			a[i][j] = a[i - 1][j - 1] + a[i - 1][j];
		}
	}
}
int main() {
	int n;
	make();
	while (scanf("%d", &n) != EOF) {
		int i, j;
		printf("1\n");
		for (i = 2; i <= n; ++i) {
			printf("1");	
			for (j = 2; j <= i; ++j) {
				printf(" %d", a[i][j]);
			}
			printf("\n");
		}
		printf("\n");
	}
	return 0;
}
