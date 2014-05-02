#include <cstdio>
#include <cstring>

char s[100000];
char t[100000];

int main(){
	while(scanf("%s %s", s, t) != EOF){
			
		int ls = strlen(s);
		int lt = strlen(t);

		int i, j, cur = -1;
		if (ls > lt)
			printf("No\n");
		else {
			int count = 0;
			for (i = 0; i < ls; ++i){
				// here we got every letter
				for (j = cur+1;j < lt; ++j){
					if (s[i] == t[j]){
						cur = j;
						count++;
						break;
					}
				}
			}
			if (count == ls)
				printf("Yes\n");
			else
				printf("No\n");
		}
	}	
	return 0;
}