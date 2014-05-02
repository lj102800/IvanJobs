#include <cstdio>
#include <cstring>

char line[1000];

int main() {
	int n;
	scanf("%d", &n);
	getchar();
	while(n--) {
		gets(line);
		int i = 0;
		int count = 0;
		int len = strlen(line);
		for( ; i < len; ++i) {
			if (line[i] < 0)
				count++;
		}	
		printf("%d\n", count / 2);	
	}

	return 0;
}
