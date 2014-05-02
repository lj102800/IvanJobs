#include <cstdio>
#include <vector>
#include <cstring>

using namespace std;

vector<int> v;
char str[15];

void printn(int n, char c) {
		int i;
		for (i = 0; i < n; ++i)
			printf("%c", c);
}
void make(char str[]) {
		v.clear();
		int n = strlen(str);
		int i;
		for (i = 0; i < n; ++i)
			v.push_back(str[i] - '0');
}

int main() {
		int s, n;
		while (scanf("%d %s", &s, str) && s!=0) {
				int rown = 2 * s + 3;
				int coln =  s + 2;
				make(str);
				int dn = v.size();
				
				int i, j;
				for (j = 1; j <= rown; ++j){
					for (i = 0; i < dn; ++i)  {
							int d = v[i];
								if (d == 0) {
									if (j == s + 2){
										printn(s + 2, ' ');
									} else if (j == 1 || j == 2 * s + 3) {
										printf(" ");
										printn(s, '-');
										printf(" ");
									} else if ((j >=2 && j <= s + 1) || (j >= s + 3 && j <= 2 * s + 2)) {
										printf("|");
										printn(s, ' ');
										printf("|");
									}
								} else if (d == 1) {
									if (j == 1 || j == s + 2 || j == 2 * s + 3) {
										printn(s + 2, ' ');
									} else if ((j >=2 && j <= s + 1) || (j >= s + 3 && j <= 2 * s + 2)) {
										printn(s + 1, ' ');
										printf("|");
									}
								} else if (d == 2) {
									if (j == 1 || j == s + 2 || j == 2 * s + 3) {
										printf(" ");
										printn(s, '-');
										printf(" ");
									} else if (j >=2 && j <= s + 1 ) {
										printn(s + 1, ' ');
										printf("|");
									} else if (j >= s + 3 && j <= 2 * s + 2) {
										printf("|");
										printn(s + 1, ' ');
									}
								} else if (d == 3) {
									if (j == 1 || j == s + 2 || j == 2 * s + 3) {
										printf(" ");
										printn(s, '-');
										printf(" ");
									} else if ((j >=2 && j <= s + 1) || (j >= s + 3 && j <= 2 * s + 2)) {
										printn(s + 1, ' ');
										printf("|");
									}
								} else if (d == 4) {
									if (j == 1 || j == 2 * s + 3) {
										printn(s + 2, ' ');
									} else if (j == s + 2) {
										printf(" ");
										printn(s, '-');
										printf(" ");
									}
									else if ((j >=2 && j <= s + 1)  ) {
										printf("|");
										printn(s, ' ');
										printf("|");
									} else if (j >= s + 3 && j <= 2 * s + 2) {
										printn(s + 1, ' ');
										printf("|");
									}
								} else if (d == 5) {
									if (j == 1 || j == s + 2 || j == 2 * s + 3) {
										printf(" ");
										printn(s , '-');
										printf(" ");
									} else if (j >= 2 && j <=s + 1) {
										printf("|");
										printn(s + 1, ' ');
									} else if (j >= s + 3 && j <= 2 * s + 2) {
										printn(s + 1, ' ');
										printf("|");
									}
								} else if (d == 6) {
									if (j == 1 || j == s + 2 || j == 2 * s + 3) {
										printf(" ");
										printn(s, '-');
										printf(" ");
									} else if ((j >= s + 3 && j <= 2 * s + 2)) {
										printf("|");
										printn(s, ' ');
										printf("|");
									} else if (j >=2 && j <= s + 1) {
										printf("|");
										printn(s + 1, ' ');
									}
								} else if (d == 7) {
									if ( j == s + 2 || j == 2 * s + 3) {
										printn(s + 2, ' ');
									} else if (j == 1) {
										printf(" ");
										printn(s , '-');
										printf(" ");
									}else if ((j >=2 && j <= s + 1) || (j >= s + 3 && j <= 2 * s + 2)) {
										printn(s + 1, ' ');
										printf("|");
									}
								} else if (d == 8) {
									if (j == 1 || j == s + 2 || j == 2 * s + 3) {
										printf(" ");
										printn(s, '-');
										printf(" ");
									} else if ((j >=2 && j <= s + 1) || (j >= s + 3 && j <= 2 * s + 2)) {
										printf("|");
										printn(s, ' ');
										printf("|");
									}
								} else if (d ==9) {
									if (j == 1 || j == s + 2 || j == 2 * s + 3) {
										printf(" ");
										printn(s, '-');
										printf(" ");
									} else if ((j >=2 && j <= s + 1) ) {
										printf("|");
										printn(s, ' ');
										printf("|");
									} else if (j >= s + 3 && j <= 2 * s + 2) {
											printn(s + 1, ' ');
											printf("|");
									}
								}
								if (i != dn - 1)
									printf(" ");
								else
									printf("\n");
					}
					
				}
				printf("\n");
		}
		return 0;
}