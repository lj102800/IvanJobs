///**
// * demo for longest common subsequence
// * one of subsequences of 1234567 is 135  
// */
//
//#include <cstdio>
//
//#define MAX_LEN 100
//#define UP 'U'
//#define LEFT 'L'
//#define UPLEFT 'V'
//
//char X[MAX_LEN];
//char Y[MAX_LEN];
//
//int c[MAX_LEN][MAX_LEN];
//char b[MAX_LEN][MAX_LEN];
//
//
//class LCS {
//public:
//	void doIt(char * fileIn, char * fileOut);
//	
//	void printLCS(char b[][MAX_LEN], char x[], int m, int n);
//};
//
//void LCS::printLCS(char b[][MAX_LEN], char x[], int m , int n) {
//	if (m == 0 || n == 0) return ;
//	
//	if (b[m][n] == UPLEFT) {
//		printLCS(b, x, m - 1, n - 1);
//		printf("%c\n", x[m]);
//	} else  if (b[m][n] == UP) {
//		printLCS(b, x, m - 1, n);
//	} else {
//		printLCS(b, x, m, n - 1);
//	}
//}
//
//void LCS::doIt(char * fileIn, char * fileOut) {
//	FILE *pFin = fopen(fileIn, "r");
//	FILE *pFout = fopen(fileOut, "w");
//
//	int n_cases;
//	fscanf(pFin, "%d\n", &n_cases);
//
//	int i;
//	for(i = 1; i <= n_cases; ++i) {
//		// in cases
//		int m, n;
//		fscanf(pFin, "%d %d\n", &m, &n);
//		int j, k;
//		for (j = 1; j <= m; ++j) {
//			if (j != m)
//				fscanf(pFin, "%c ", &X[j]);
//			else 
//				fscanf(pFin, "%c\n", &X[j]);
//		}
//		for (k = 1;k <= n; ++k) {
//			if (k != n)
//				fscanf(pFin, "%c ", &Y[k]);
//			else
//				fscanf(pFin, "%c\n", &Y[k]);
//		}
//
//		//make zeros
//		int i;
//		for ( i = 0; i <= m; ++i ) c[i][0] = 0;
//		for ( i = 1; i <= n; ++i) c[0][i] = 0;
//
//		//dynamic programming with iterations
//		for (i = 1; i <= m; ++i) {
//			for (j = 1; j <= n; ++j) {
//				if (X[i] == Y[j]) {
//					c[i][j] = c[i - 1][j - 1] + 1;
//					b[i][j] = UPLEFT;
//				} else if (c[i - 1][j] >= c[i][j - 1]) {
//					c[i][j] = c[i - 1][j];
//					b[i][j] = UP;
//				} else {
//					c[i][j] = c[i][j - 1];
//					b[i][j] = LEFT;
//				}
//			}
//		}
//		printLCS(b, X, m , n);
//		//before ending cases' loop
//	}
//	fclose(pFin);
//	fclose(pFout);
//}
//
////int main()
////{
////    LCS i;
////    i.doIt("in.txt", "out.txt");
////
////    return 0;
////}