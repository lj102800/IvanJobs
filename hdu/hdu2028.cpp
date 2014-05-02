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

int lcm(int a, int b) {
	int gcd = euclid_gcd(a, b);
	int res = gcd;
	while (gcd != 1) {
		a /= gcd;
		b /= gcd;
		gcd = euclid_gcd(a, b);		
	}
	res *= (a * b);
	return res;
}

int main() {
	int n;
	int curr, val;
	while (scanf("%d", &n) != EOF) {
		scanf("%d", &val);
		n--;
		curr = val;
		while (n--) {
			scanf("%d", &val);
			curr = 	lcm(curr, val);	
		}
		printf("%d\n", curr);
	}
	return 0;
}
