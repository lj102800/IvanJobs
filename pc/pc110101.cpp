#include <cstdio>
#include <map>

using namespace std;

map<int, int> m; // used for mapping n and its cycle length.
map<int, int>::iterator it;

int  make(long long i) {
	if (i == 1)
		return 1;
	it = m.find(i);
	if (it != m.end()) {
			return it->second;
	} else {
			long long oldi = i;
			if (i % 2 == 0) i /= 2;
			else i = 3 * i + 1;
			int v = 1 + make(i);
			m.insert(pair<int,  int>(oldi, v));
			return v;
	}
}

int main() {

	int a, b;
	while (scanf("%d%d", &a, &b) != EOF) {
		
		int c, d;
		c = a;
		d = b;
		if (a > b) {
				int tmp = a;
				a = b;
				b = tmp;
		}
		
		long long i;
		int max = -1;
		for (i = a; i <= b; ++i) {
				int len = make(i);
				if (len > max) max = len;
		}
		
		printf("%d %d %d\n", c, d, max);
		
	}
	return 0;
}
