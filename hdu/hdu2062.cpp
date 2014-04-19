#include <iostream>
#include <cstring>
#include <vector>
#include <cstdio>

using namespace std;
#define L(fmt, ...) do {if (false) printf(fmt "\n", ##__VA_ARGS__);}while(false)
long long N[25];
int V[25];
vector<int> vi;

void makeN() {
    N[1] = 1;
    int i;
    for (i = 2; i <= 20; ++i)
        N[i] = i * (1 + N[i - 1]);
}

int seekMin(int n) {
    int i;
    int cnt = 0;
    for (i = 1; i < 25; ++i) {
        if (V[i] == 0) {
                cnt++;
                if (cnt == n){
                    L("seekMin(%d)=%d", n, i);
                    V[i] = 1;
                    return i;
                }
        }
    }
    return -1;
}

void make(int n, long long m) {
    L("make(%d, %lld)", n, m);
    if (m < 0 ) return ;
    if (n == 1 && m == 0) {
        vi.push_back(seekMin(1));
        return ;
    }
    long long k = (1 + N[n - 1]);
    long long a, b;
    a = m / k; b = m % k;
    vi.push_back(seekMin((int)(a + 1)));
    make(n - 1, b - 1);
}

void display() {
    int len = vi.size();
    int i;
    cout<<vi[0];
    for (i = 1; i < len; ++i)
        cout<<" "<<vi[i];
    cout<<endl;
}

int main() {
    ios::sync_with_stdio(false);
    int n;
    long long m;
    makeN();
    while(cin>>n>>m) {
        vi.clear();
        memset(V, 0, sizeof(V));
        make(n, m - 1);
        display();
    }
    return 0;
}

