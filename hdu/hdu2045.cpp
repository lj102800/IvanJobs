#include <cstdio>

long long s[60];

void build() {
	s[1] = 3;
	s[2] = 6;
	s[3] = 6;
	int i;
	for (i = 4; i < 60; ++i) {
		s[i] = s[i - 1] + 2 * s[i - 2];
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
