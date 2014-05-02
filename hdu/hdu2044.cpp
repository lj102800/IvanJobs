#include <cstdio>

long long fb[70];

void build() {
	fb[1] = 1;
	fb[2] = 1;
	int i;
	for (i = 3; i < 70; ++i)
		fb[i] = fb[i - 1] + fb[i - 2];
}

int main() {
	build();
	int n, a, b;
	scanf("%d", &n);
	while (n--){
		scanf("%d%d", &a, &b);
		printf("%I64d\n", fb[b - a + 1]);
	}
	return 0;
}
