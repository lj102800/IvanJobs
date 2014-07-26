#include <cstdio>
typedef long long LL;
LL euclid_gcd(LL m, LL n){
	LL r = m % n;
	while(r != 0){
		m = n;
		n = r;
		r = m % n;
	}
	return n;
}

// find int x, y, make ax + by = d, d = gcd(a, b) and Minimize |x| + |y|.
// although a and b are int, x and y may bigger than int
void gcd(LL a, LL b, LL& d, LL& x, LL& y) {
	if (!b) { d = a; x = 1; y = 0;}
	else{ gcd(b, a % b, d, y , x); y -= x * (a / b); }
}

int main() {
	LL m, n;
	while (scanf("%lld%lld", &m, &n) && !(m == 0 && n == 0)){
		printf("%lld\n", euclid_gcd(m, n));
	}

	return 0;
}
