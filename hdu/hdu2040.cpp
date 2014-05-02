#include <cstdio>

int main() {
	int m, a, b;
	scanf("%d", &m);
	while (m--) {
		scanf("%d%d", &a, &b);
		int i;
		int s1, s2;
		s1 = s2 = 0;
		for (i = 1; i < a; ++i) {
			if (a % i == 0)
				s1 += i;
		}
		for (i = 1; i < b; ++i) {
			if (b % i == 0) 
				s2 += i;
		}
		if (s1 == b && s2 == a)
			printf("YES\n");
		else
			printf("NO\n");
	}
	return 0;
}
