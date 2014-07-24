/*
	T: text string
	P: pattern string
	f: fail array
 */
void kmp(char* T, char* P, int* f) {
	int n = strlen(T), m = strlen(P);
	getFail(P, f);
	int j = 0;
	for (i = 0; i < n; i++) {
		while(j && P[j] != T[j]) j = f[j];
		if(P[j] == T[j]) j++;
		if (j == m) printf("%d\n", i - m + 1); // find one match here
	}
}

void getFail(char* P, int* f) {
	int m = strlen(P);
	f[0] = 0; f[1] = 0;
	for (int i = 1; i < m; i++) {
		int j = f[i];
		while(j && P[i] != P[j]) j = f[j];
		f[i + 1] = P[i] == P[j] ? j + 1 : 0;
	}
}
