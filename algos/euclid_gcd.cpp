#include <cstdio>

int euclid_gcd(int m, int n){
	int r = m % n;
	while(r != 0){
		m = n;
		n = r;
		r = m % n;
	}
	return n;
}

int main() {
	int m, n;
	while (scanf("%d%d", &m, &n) && !(m == 0 && n == 0)){
		printf("%d\n", euclid_gcd(m, n));
	}

	return 0;
}