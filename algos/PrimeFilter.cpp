const int MAXN = 10000000 + 10
const int MAXP = 700000;

int vis[MAXN];
int prime[MAXP];

void sieve(int n) {
	int m = (int) sqrt(n + 0.5);
	memset(vis, false, sizeof(vis));
	for (int i = 2; i <= m; i++) if(!vis[i]) {
		for(j = i * i; j <= n; j += i) vis[j] = true;
	}
}

int genPrimes(int n) {
	sieve(n);
	int c = 0;
	for(int i = 2; i <= n; i++) if (!vis[i]) prime[c++] = i;
	return c;
}
