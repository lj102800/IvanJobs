#include <cstdio>
#include <cstring>

char pass[100];
char ch[] = {'~', '!', '@', '#', '$', '%', '^'};
bool isin(char c) {
	int i;
	for (i = 0; i < sizeof(ch); ++i) {
		if (c == ch[i])
			return true;
	}
	return false;
}

int main() {
	int m;
	scanf("%d", &m);
	while (m--) {
		scanf("%s", pass);
		bool b1, b2, b3, b4;
		b1 = b2 = b3 = b4 = false;
		int cnt = 0;
		int len = strlen(pass);
		int i;
		for (i = 0; i < len; ++i) {
			if (b1 == false && pass[i] >= 'A' && pass[i] <= 'Z') {
				b1 = true;
				cnt++;
			}
			if (b2 == false && pass[i] >= 'a' && pass[i] <= 'z') {
				b2 = true;
				cnt++;
			}
			if (b3 == false && pass[i] >= '0' && pass[i] <= '9') {
				b3 = true;
				cnt++;
			}
			if (b4 == false && isin(pass[i])) {
				b4 = true;
				cnt++;
			}
		} 
		if (cnt >= 3 && len >= 8 && len <= 16)
			printf("YES\n");
		else
			printf("NO\n");
	}
	return 0;
}
