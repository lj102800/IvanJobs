#include <cstdio>
#include <cctype>
#include <cstring>

char line[110];
int main() {
	while (gets(line)) {
		int i;
		bool blank = true;
		for (i = 0; i < strlen(line); ++i) {
			if (blank && line[i] != ' ') {
				line[i] = toupper(line[i]);
				blank = false;
			} 
			if (!blank && line[i] == ' ') {
				blank = true;
			}
			
		}
		printf("%s\n", line);
	}	
	return 0;
}
