#include <cstdio>
#include <cstring>

#define MAX 200
char mat[MAX][MAX];
int rec[MAX][MAX];
int d[][2] = {{-1, -1}, {-1, 0}, {-1, 1}, {0, 1}, {1, 1}, {1, 0}, {1, -1}, {0, -1}};

void make(int n, int m) {
	int i, j;
	for (i = 0; i < n; ++i) {
			for (j = 0; j < m; ++j) {
					if (mat[i][j] == '*') {
							rec[i][j] = -1;
							int k;
							for (k = 0; k < 8; ++k) {
									int newi = i + d[k][0];
									int newj = j + d[k][1];
									if (newi >= 0 && newi < n && newj >= 0 && newj < m) {
										if (mat[newi][newj] != '*')
											rec[newi][newj]++;
									}
							}
					}
			}
	}
}

void display(int n, int m, int fieldid) {
	if (fieldid > 1)
		printf("\n");
	printf("Field #%d:\n", fieldid);
	int i, j;
	for (i = 0; i < n; ++i) {
			for (j = 0; j < m; ++j) {
					if (rec[i][j] != -1)
						printf("%d", rec[i][j]);
					else
						printf("*");
			}
			printf("\n");
	}

}

int main(int argc, char **argv)
{
	int n, m;
	scanf("%d %d", &n, &m);
	getchar();
	int fieldid = 0;
	while (!(n == 0 && m == 0)) {
			memset(rec, 0, sizeof rec);
			int i, j;
			for (i = 0; i < n; ++i) {
				for (j = 0; j < m; ++j)
					scanf("%c", &mat[i][j]);
				getchar();
			}
			
			make(n, m);
			display(n, m, ++fieldid); // don't use fieldid++ here, it's wrong.
			
		scanf("%d %d", &n, &m);
		getchar();
	}
	return 0;
}
