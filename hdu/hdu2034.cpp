#include <iostream>
#include <list>
#include <algorithm>
using namespace std;
void p(int i) {
	cout<<i<<" ";
}
list<int> a, b;
int main() {
	int n, m;
	while (cin>>n>>m && !(n == 0 && m == 0)) {
		int i, val;
		a.clear();
		b.clear();
		for (i = 0; i < n; ++i) {
			cin>>val;
			a.push_back(val);
		}
		for (i = 0; i < m; ++i) {
			cin>>val;
			b.push_back(val);
		}
		
		list<int>::iterator it1, it2;
		it1 = a.begin();
		while (it1 != a.end()) {
			it2 = find(b.begin(), b.end(), *it1);
			if (it2 != b.end()){
				it1 = a.erase(it1);
			} else {
				it1++;
			}		
		}

		if (a.empty())
			cout<<"NULL"<<endl;
		else {
			a.sort();
			for_each(a.begin(), a.end(), p);
			cout<<endl;	
		}
	}
	return 0;
}
