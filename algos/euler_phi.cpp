// phi(n), number of integers not bigger than n and prime with n
int eulerPhi(int n) {
	int m = (int)sqrt(n + 0.5);
	int ans = n;
	for(int i = 2; i <= m; i++) if (n % i == 0) {
		ans = ans / i * (i - 1);
		while(n % i == 0) n /= i;
	}
	if (n > 1) ans = ans / n * (n - 1);
}

// similar with primers' sieve, to calculate phi(1), phi(2), phi(3)...
#define MAXN 1000
int phi[MAXN];
void phiTable(int n) {
	for(int i = 2; i <= n; i++) phi[i] = 0;
	phi[1] = 1;
	for(int i = 2; i <= n; i++) if(!phi[i])
		for(int j = i; j <= n; j += i) {
			if(!phi[j]) phi[j] = j;
			phi[j] = phi[j] / i * (i - 1);
		}
}
