#include <cstdio>
#include <cstring>
#include <cmath>

char line[30];
char a[10];
char b[10];
char c[10];

int main(){
	while(scanf("%s", line)){
		
		char *pchr = NULL;
		int i, m, n, len;
		len = strlen(line);
		pchr = strchr(line, '+');
		m = pchr - line;
		pchr = strchr(line, '=');
		n = pchr - line;

		int suma = 0;
		for (i = m-1; i >= 0; i--){
			suma += (line[i] - '0') * (int)pow(10.0, i * 1.0);
		}

		int sumb = 0;
		for (i = n-1; i >= m+1; i--){
			sumb += (line[i] - '0') * (int)pow(10.0, i * 1.0 - m - 1);
		}

		int sumc = 0;
		for (i = len -1; i >= n+1; i--){
			sumc += (line[i] - '0') * (int)pow(10.0, i * 1.0 - n - 1);
		}

		if (suma + sumb == sumc){
			printf("True\n");
		}else{
			printf("False\n");
		}
		if(strcmp(line, "0+0=0") == 0)
			break;
		
	}
	return 0;
}