#include <cstdio>

int main() {
	int a, b;
	while (scanf("%d%d", &a, &b) && !(a == 0 && b == 0)) {
		int res = 1;
		int i;
		for (i = 1; i <= b; ++i) {	
			res = ((res % 1000) * (a % 1000)) % 1000;	
		}
		printf("%d\n", res);
	}
	return 0;
}
