#include <cstdio>
int fb[50];
void make() {
	fb[1] = 1;
	fb[2] = 1;
	int i;
	for (i = 3; i < 50; ++i)
		fb[i] = fb[i - 1] + fb[i - 2];
}
int main() {
	make();
	int n;
	scanf("%d", &n);
	while (n--) {
		int val;
		scanf("%d", &val);
		printf("%d\n", fb[val]);

	}
	return 0;
}
