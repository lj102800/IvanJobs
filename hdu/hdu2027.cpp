#include <cstdio>
#include <cstring>

char line[110];

int main() {
	int n;
	scanf("%d", &n);
	int m = n;
	getchar();
	while (n--) {
		gets(line);
		int len = strlen(line);
		int i;
		int a, e, ii, o, u;
		a = e = ii = o = u = 0;
		for (i = 0; i < len; ++i) {
			switch(line[i]) {
				case 'a':
					a++;
					break;
				case 'e':
					e++;
					break;
				case 'i':
					ii++;
					break;
				case 'o':
					o++;
					break;
				case 'u':
					u++;
					break;
				default:
					break;
			}
		}
		if (n != m - 1)
			printf("\n");
		printf("a:%d\ne:%d\ni:%d\no:%d\nu:%d\n", a, e, ii, o, u);
				
	}
	return 0;
}
