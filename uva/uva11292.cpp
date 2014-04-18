#include <iostream>
#include <algorithm>
using namespace std;
int H[20010];
int P[20010];

int main() {
	ios::sync_with_stdio(false);
	int n, m, i, j;
	while (cin>>n>>m && !(n == 0 && m == 0)) {
		for (i = 0; i < n; ++i)
			cin>>H[i];
		for (i = 0; i < m; ++i)
			cin>>P[i];
		sort(H, H + n);
		sort(P, P + m);
		i = j = 0;
		int res = 0;
		do {
			if(i == n || j == m) break;
			if (P[j] >= H[i]) {
				res += P[j];
				i++; j++;
			} else {
				j++;
			}
		} while(true);
		if (i != n) cout<<"Loowater is doomed!"<<endl;
		else cout<<res<<endl;	
	}
	return 0;
}
