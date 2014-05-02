#include <cstdio>

int gap[] = {100, 50, 10, 5, 2, 1};
int make(int val) {
	int count = 0;
	int i;
	for (i = 0; i < 6; ++i) {
		if (val >= gap[i]) {
			count += val / gap[i];
			val -= gap[i] * (val / gap[i]);	
		}
	}
	return count;
}
int main() {
	int n;
	while (scanf("%d", &n) && n != 0) {
		int i;
		int count = 0;
		int val;
		for (i = 0; i < n; ++i) {
			scanf("%d", &val);
			count += make(val);
		}
		printf("%d\n", count);	
	}
	return 0;
}
