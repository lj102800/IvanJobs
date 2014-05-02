#include <cstdio>

int main() {
	int n;
	int ah, am, as, bh, bm, bs;
	scanf("%d", &n);
	while(n--) {
		scanf("%d%d%d%d%d%d", &ah, &am, &as, &bh, &bm, &bs);
		as += bs;
		am += (as / 60);
		as %= 60;
		am += bm;
		ah += (am / 60);
		am %= 60;
		ah += bh;
		printf("%d %d %d\n", ah, am, as);
	}
	return 0;
}
