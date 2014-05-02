#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
	int n, m;
	while (cin>>n>>m && !(n == 0 && m == 0))  {
		int i;
		vector<int> v;
		v.clear();
		v.push_back(m);
		int tmp;
		for (i = 0; i < n; ++i) {
			cin>>tmp;
			v.push_back(tmp);
		}
		sort(v.begin(), v.end());
		cout<<v[0];
		for (i = 1; i < v.size(); ++i)
			cout<<" "<<v[i];
		cout<<endl;
	}
	return 0;
}
