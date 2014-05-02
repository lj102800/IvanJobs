#include <cstdio>
struct Lpoint{ double x, y;};

struct Lpoint v[110];

double areaofp(int vcount, Lpoint plg[]) {
	int i;
	double s;
	if (vcount < 3) {
		return 0;
	}

	s = plg[0].y * (plg[vcount - 1].x - plg[1].x);
	for (i = 1; i < vcount; i++) {
		s += plg[i].y * (plg[i - 1].x - plg[(i + 1) % vcount].x);
	}
	return s/2;
}
int main() {
	int n;
	while (scanf("%d", &n) && n != 0) {
		int i;
		for (i = 0; i < n; i+=1) {
			scanf("%lf", &(v[i].x));
			scanf("%lf", &(v[i].y));
		}

		printf("%.1lf\n", areaofp(n, v));
	}
	return 0;
}
