#include <cstdio>

#define MAX 10010
float mat[MAX];

float max_float (int n) {
	int i;
	float res = -1.0;
	for (i = 0; i < n; ++i) {
		if (mat[i] > res)
			res = mat[i];
	}
	return res;
}

int num (float c, int n) {
	int res = 0;
	int i;
	for (i = 0; i < n; ++i) {
		res += int(mat[i] / c);
	}
	return res;
} 

int main() {
	int n, N, K;
	int sh, i;
	float res;
	scanf("%d", &n);
	for (sh = 0; sh < n; ++sh){
		scanf("%d%d", &N, &K);
		for (i = 0; i < N; ++i) {
			scanf("%f", &mat[i]);
		}
		if (sh != 0)
			printf("\n");
		float up_end = max_float(N);
		float low_end = 0.0;
		float mid = (up_end + low_end) / 2;
		float res;
		int curr_n = num(up_end, N);
		if (curr_n >= K) {
			printf("%.2f\n", up_end - 0.005);
		} else {
			do {
				curr_n = num(mid, N);
				if (curr_n > K) {
					low_end = mid;	
				} 
				if (curr_n == K){
					res = mid;
				}
				if (curr_n < K) {
					up_end = mid;
				}
				mid = (up_end + low_end) / 2;
			} while (curr_n != K);
		}
		
		
		if (res >= 0.01) {
			printf("%.2f\n", res-0.005);
		} else {
			printf("0.00\n");
		}
	}

	return 0;
}