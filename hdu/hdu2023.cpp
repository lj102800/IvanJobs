#include <cstdio>

double mat[60][10];
double averStu[60];
double averSub[10];

int main() {
	int n, m, i, j, count;
	double su;
	while (scanf("%d%d", &n, &m) != EOF) {
		for (i = 0; i < n; ++i) {
			su = 0;
			for (j = 0; j < m; ++j) {
				scanf("%lf", &(mat[i][j]));
				su += mat[i][j];
			}
			averStu[i] = su / m;		
		}
		
		for (i = 0; i < m; ++i) {
			su = 0;
			for (j = 0; j < n; ++j) {
				su += mat[j][i];	
			}
			averSub[i] = su / n;
		}

		count = 0;
		bool isc = true;
		for (i = 0; i < n; ++i) {
			isc = true;
			for (j = 0; j < m; ++j) {
				if (mat[i][j] < averSub[j]) {
					isc = false;
					break;
				}
					
			}
			if (isc) {
				count++;
			}
		}
		printf("%.2lf", averStu[0]);
		for (i = 1; i < n; ++i)
			printf(" %.2lf", averStu[i]);
		printf("\n");
		printf("%.2lf", averSub[0]);
		for (i = 1; i < m; ++i)
			printf(" %.2lf", averSub[i]);
		printf("\n");
		printf("%d\n\n", count); 
			
		
	}	
	return 0;
}
