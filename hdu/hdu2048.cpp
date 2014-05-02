#include <cstdio>
double s[30];
void build() {
	s[1] = 0;
	s[2] = 1;
	int i;
	for (i = 3; i < 30; ++i) {
		s[i] = (i - 1) * (s[i - 1] + s[i - 2]);
	}
}

double fac(int n) {
	double res = 1;
	int i;
	for (i = 1; i <=n; ++i)
		res *= i;
	return res;
}
int main() {
	build();
	int C;
	scanf("%d", &C);
	while (C--) {
		int n;
		printf("%.2lf%%\n", (double)s[n] * 100 / (double)fac(n));
	}
	return 0;
}
