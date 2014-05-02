#include <iostream>
#include <list>

struct TSeg {int s, e;};
using namespace std;

list<TSeg> li;

bool cmp(TSeg t1, TSeg t2) {
	if (t1.s < t2.s)
		return true;
	else if (t1.s > t2.s) 
		return false;
	else {
		if (t1.e < t2.e)
			return true;
		else {
			return false;
		}
	}	
}

int main() {
	int n;
	while (cin>>n && n != 0) {
		int i;
		li.clear();
		for (i = 0; i < n; ++i) {
			TSeg t;
			cin>>t.s>>t.e;
			li.push_back(t);
		}
		li.sort(cmp);

		list<TSeg>::iterator it1, it2;
		it1 = li.begin();
		it2 = li.begin();
		it2++;
		while (it2 != li.end()) {
			TSeg g1 = *it1;
			TSeg g2 = *it2;
			if (g2.e <= g1.e) {
				it1 = li.erase(it1);
				it2++;
			} else if (g2.s >= g1.e) {
				it1++;
				it2++;
			} else {
				it2 = li.erase(it2);
			}
		}
		cout<<li.size()<<endl;		
	}
	return 0;
}
