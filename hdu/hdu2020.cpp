#include <iostream>
#include <algorithm>
#include <vector>
#include <cstdlib>
using namespace std;

bool AbsCmp(int a, int b) {
	return abs(b) < abs(a);
}

vector<int> v;
int main() {
	int n;
	while (cin>>n && n != 0) {
		v.clear();
		int tmp;
		int i;
		for (i = 0; i < n; ++i) {
			cin>>tmp;
			v.push_back(tmp);
		}
		sort(v.begin(), v.end(), AbsCmp);
		cout<<v[0];
		for (i = 1; i < n; ++i) {
			cout<<" "<<v[i];
		}	
		cout<<endl;
	}
	return 0;
}
