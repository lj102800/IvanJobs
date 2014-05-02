#include <cstdio>

int main() {
	int n;
	float a, b, c;
	scanf("%d", &n);
	while (n--) {
		scanf("%f%f%f", &a, &b, &c);
		bool itcan = true;
		if (a + b <= c)
			itcan = false;
		if (a + c <= b)
			itcan = false;
		if (b + c <= a)
			itcan = false;
		if (itcan)
			printf("YES\n");
		else 
			printf("NO\n");
				
	}
	return 0;
}
