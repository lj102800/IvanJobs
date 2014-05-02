#include <cstdio>
#include <vector>

using namespace std;

vector<int> v;

int main() {
	int n, m;
	
	// 0, 11, 12,  21, 22
	int s = 0;
	

	while(scanf("%d%d", &n, &m) != EOF) {
		int i = 1;
		int count = 0;
		int su = 0;
		v.clear();
		while ( i <= n ) {
			su += i * 2;
			count++;
			if (count == m) {
				v.push_back(su / m);
				su = 0;
				count = 0;			
			}
			i++;	
		}
		if (count < m && count > 0){
			v.push_back(su / count);
		}

		printf("%d", v[0]);
		for (i = 1; i < v.size(); ++i) {
			printf(" %d", v[i]);
		}
		printf("\n");
	
	}
	return 0;
}
