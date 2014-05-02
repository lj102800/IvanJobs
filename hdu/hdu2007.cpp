//#include <cstdio>
//
//#define C(i) i*i*i
//#define S(i) i*i
//int main() {
//	int m, n;
//	while(scanf("%d%d", &m, &n) == 2) {
//		int i;
//		int resSquare = 0;
//		int resCubic = 0;
//		if (m > n) {
//			int tmp = m;
//			m = n;
//			n = tmp;
//		}
//		for (i = m; i <= n; ++i) {
//			if (i % 2 == 1)
//				resCubic += C(i);
//			else
//				resSquare += S(i);
//		}
//
//		printf("%d %d\n", resSquare, resCubic);
//	}
//
//	return 0;
//}