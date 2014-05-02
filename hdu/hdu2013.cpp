#include <cstdio>

int make(int n) {
	int s = 1;
	while (n != 1) {
		s = (s + 1) * 2;
		n--;		
	}
	return s;

}
int main() {
	int n;
	while (scanf("%d", &n) != EOF) {
		printf("%d\n", make(n));
	}
	return 0;
}
