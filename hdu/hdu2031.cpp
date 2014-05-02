#include <cstdio>
#include <vector>

using namespace std;


char letter[] = {'0', '1', '2', '3', '4',
		 '5', '6', '7', '8', '9',
		 'A', 'B', 'C', 'D', 'E',
		 'F'};
int main() {
	int N, R;
	bool negative = false;
	while (scanf("%d%d", &N, &R) != EOF) {
		negative = false;
		if (N < 0) {
			negative = true;
			N = -N;
		}

		if (negative)
			printf("-");
		vector<char> v;
		v.clear();
		while (N != 0) {
			v.push_back(letter[N % R]);
			N /= R;
		}
		int i;
		for (i = v.size() - 1; i >= 0; i--) {
			printf("%c", v[i]);
		}
		printf("\n");		
	} 
	return 0;
}
