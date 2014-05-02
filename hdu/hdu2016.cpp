#include <cstdio>

int mat[130];

int main() {
	int n;
	while(scanf("%d", &n) && n != 0) {
		int i;
		for (i = 0; i < n; ++i) {
			scanf("%d", &(mat[i]));	
		}

		int minIdx = 0;
		int minVal = mat[0];
		for (i = 1; i < n; ++i) {
			if (minVal > mat[i]) {
				minVal = mat[i];
				minIdx = i;
			}
		}	
		
		int tmp = mat[minIdx];
		mat[minIdx] = mat[0];
		mat[0] = tmp;
		printf("%d", mat[0]);
		for (i = 1; i < n; ++i) {
			printf(" %d", mat[i]);
		}
		printf("\n");
	}
	return 0;
}
