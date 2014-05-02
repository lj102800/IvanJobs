//#include <cstdio>
//
//bool isPrime(int y) {
//	if (y % 4 == 0 && y % 100 != 0)
//		return true;
//	if (y % 400 == 0)
//		return true;
//	return false;
//}
//
//int month[] = {	0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
//
//int main() {
//	int y, m, d;
//
//	while (scanf("%d/%d/%d", &y, &m, &d) == 3) {
//		if (isPrime(y)) month[2] = 29;
//		else month[2] = 28;
//
//		int i, j;
//		int counter = 0;
//		for (i = 1; i <m; ++i) {
//			counter += month[i];
//		}
//		counter += d;
//		printf("%d\n", counter);
//	}
//	return 0;
//}