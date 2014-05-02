#include <cstdio>
#include <cmath>

#define PI 3.1415926

int make(float x, float y) {
	float r = sqrt(x*x + y*y);
	float area = 0.5*PI * r * r;
	return (int)area/50 + 1;
}

int main(){
	int n;
	scanf("%d", &n);
	float x, y;
	int i;
	for (i = 0; i < n; ++i) {
		scanf("%f%f", &x, &y);
		printf("Property %d: This property will begin eroding in year %d.\n", i + 1, make(x, y));
	}
	printf("END OF OUTPUT.\n");
	return 0;
}
