#include <cstdio>

long long s[60];
void build() {
	s[1] = 1;
	s[2] = 2;
	int i;
	for (i = 3; i < 60; ++i)
		s[i] = s[i - 1] + s[i - 2];
}

int main() {
	build();
	int n;
	while (scanf("%d", &n) != EOF) {
		printf("%I64d\n", s[n]);
	}
	return 0;
}
