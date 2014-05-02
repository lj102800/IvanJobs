#include <cstdio>
long long s[43];

void build() {
	s[1] = 3;
	s[2] = 8;
	int i;
	for (i = 3; i < 43; ++i) {
		s[i] = 2 * (s[i - 1] + s[i - 2]);
	}
}
int main() {
	build();
	int n;
	while (scanf("%d", &n) != EOF) {
		printf("%I64d\n", s[n]);
	}
	return 0;
}
