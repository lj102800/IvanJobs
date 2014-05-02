#include <cstdio>

int main() {
	int n;
	scanf("%d", &n);
	while (n--) {
		int val;
		scanf("%d", &val);
		int i;
		int res = 3;
		while (val--) {
			res = (res - 1) * 2;	
		}
		printf("%d\n", res);
	}
	return 0;
}
