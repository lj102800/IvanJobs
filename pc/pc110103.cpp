#include <cstdio>

#define MAX 1010

double cost[MAX];

double make(int n, double aver) {
	int i;
	double res1 = 0.00;
	double res2 = 0.00;
	for (i = 0; i < n; ++i) {
			if (cost[i] < aver) {
				res1 += aver - cost[i];
			} else {
				res2 += cost[i] - aver;
			}
	}
	
	return res1 < res2 ? res1 : res2;
}

int main() {
		int n;
		while (scanf("%d", &n) && n != 0) {
			int i;
			double res = 0.0;
			for (i = 0; i < n; ++i) {
					scanf("%lf", &cost[i] );
					res += cost[i];
			}
			
			double aver = res / n;
			aver = (int)(aver * 100 + 0.5) / (double) 100;
			printf("$%.2lf\n", make(n, aver));
		}
		return 0;
}