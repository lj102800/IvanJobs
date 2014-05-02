#include <cstdio>

int main() {
	int i;
	double sum = 0.00;
	double temp;
	for (i = 0; i < 12; ++i) {
		scanf("%lf", &temp);
		sum += temp;
	}

	printf("$%.2lf", sum/12);
	return 0;
}
