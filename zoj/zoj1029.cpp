#include <cstdio>
#include <list>
#include <algorithm>

using namespace std;



struct Segment {
	int x, y;
};

bool cmp( Segment s1, Segment s2) {
	return s1.x <= s2.x;
}

list<Segment> L;



int main() {
	int T, tc, N, i, x, y;
	scanf("%d", &T);
	for (tc = 1; tc <= T; ++tc) {
		L.clear();
		scanf("%d", &N);
		for (i = 0; i < N; ++i) {
			Segment seg;
			scanf("%d%d", &(seg.x), &(seg.y));
			L.push_back(seg);
		}
		
		L.sort(cmp);	
		int cnt = 0;
		while (!L.empty()) {
			cnt++;
			list<Segment>::iterator it = L.begin();	// get the fist element
			Segment curr;
			curr.x = it->x;
			curr.y = it->y;	// copy the first element to curr
			
			it = L.erase(it); // erase it from the list
			for (; it != L.end(); ) {
				if (it->x < curr.y && it->y < curr.y) {
					x = curr.x;
					y = curr.y;
					curr.x = it->x;
					curr.y = it->y;
					it->x = x;
					it->y = y;
					it++;
				} else if (it->x > curr.y) {
					curr.x = it->x;
					curr.y = it->y;
					it = L.erase(it);
				} else {
					it++;
				}
			}
			L.sort(cmp);
		}
		// here we get cnt
		printf("%d\n", 10 * cnt);
	}

	return 0;
}