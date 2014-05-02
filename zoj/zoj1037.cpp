/*
 * zoj1037.cpp
 *
 *  Created on: 2013Äê7ÔÂ3ÈÕ
 *      Author: Administrator
 */

#include <cstdio>

int main() {
	int c;
	scanf("%d", &c);
	int i, m, n;
	for(i = 1; i <= c; ++i) {
		scanf("%d%d", &m, &n);
		float res = 0.00;
		if (m % 2 ==0 || n % 2 == 0)
			res = m * n;
		else {
			res = m * n + 0.41;
		}

		printf("Scenario #%d:\n", i);
		printf("%.2f\n\n", res);
	}
	return 0;
}



