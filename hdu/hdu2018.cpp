#include <cstdio>

int mat[60];
void make() {
	mat[0] = 1;
	mat[1] = 2;
	mat[2] = 3;
	mat[3] = 4;
	int i;
	for (i = 4; i < 60; ++i) 
		mat[i] = mat[i - 1] + mat[i - 3];
}
int main() {
	int n;
	make();
	while(scanf("%d", &n) && n != 0) {
		printf("%d\n", mat[n - 1]);
	}
	return 0;
}
