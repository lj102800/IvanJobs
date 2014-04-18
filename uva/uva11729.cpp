#include <iostream>
#include <algorithm>

using namespace std;

typedef struct _node {
	_node():b(-1), j(-1) {}
	_node(int b_, int j_ ):b(b_), j(j_) {
	}
	int b;
	int j;
} node;

struct _cmp {
	/*
	Binary function that accepts two elements, return bool.
	It means:
	whether the first ele isconsidered to go before
	the second?
	*/
	bool operator() (node a, node b) {
		if (a.j > b.j) return true;
		return false; 
	}
} cmp;

node X[10000 + 10];
int main() {
	ios::sync_with_stdio(false);
	int n, b, j, i, res, bs, tn;
	tn = 0;
	while(cin>>n && n != 0) {
		tn++;
		for (i = 0; i < n; ++i)
			cin>>X[i].b>>X[i].j;
		
		sort(X, X+n, cmp);
		bs = res = 0;
		for (i = 0; i < n; ++i) {
			int newend = bs + X[i].b + X[i].j;
			res = (newend > res) ? newend : res;
			bs += X[i].b;	
		}
		cout<<"Case "<<tn<<":"<<" "<<res<<endl;		
	}
	return 0;
}
