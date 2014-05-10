#include <cstring>
#include <cmath>
#include <cstdio>

bool is[10100]; int prm[10100];

int getprm(int n) {
	int i, j, k = 0;
	int s, e = (int) (sqrt(0.0 + n) + 1);
	memset(is, 1, sizeof(is));
	prm[k++] = 2; is[0] = is[1] = 0;
	for (i = 4; i < n; i += 2) is[i] = 0;
	for (i = 3; i < e; i += 2) if (is[i]) {
		prm[k++] = i;
		for (s = i * 2, j = i * i; j < n; j += s) is[j] = 0;
	}
	for (; i < n; i += 2) if (is[i]) prm[k++] = i;
	return k;
}

int main() {
	int n = getprm(100);
	printf("%d\n", n);
	for (int i = 0; i < n; ++i)
		printf("%d\n", prm[i]);
	return 0;
}
