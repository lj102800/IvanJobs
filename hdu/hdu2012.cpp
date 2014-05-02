#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <list>
#include <cstring>

using namespace std;
list<int> s;
char mat[3000];

void make() {
	memset(mat, 0, sizeof(mat));
	int i; 
	for (i = 2; i < 3000; ++i) {
		s.push_back(i);
	}
	while (!s.empty()) {
		list<int>::iterator it = min_element(s.begin(), s.end());
		int candy = *it;
		mat[candy] = 1;
		it = s.begin();
		while (it != s.end()) {
			int item = *it;
			if (item % candy == 0) {
				it = s.erase(it);
			} else {
				it++;
			}
		}
	}
					
}

int main() {
	int x, y;
	make();
	while(scanf("%d%d", &x, &y) && !(x == 0 && y == 0)){
		x = abs(x);
		y = abs(y);
		int start = min(x, y);
		int end = max(x, y);
		int i, val;
		bool isok = true;
		for (i = start; i <= end; ++i) {
			val = i * i + i + 41;
			if (mat[val] == 0) {
				printf("Sorry\n");
				isok = false;
				break;
			}		
		}
		if (isok) {
			printf("OK\n");
		} 		
	}
	return 0;
}
