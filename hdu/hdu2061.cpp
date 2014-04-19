#include <iostream>
#include <string>
#include <iomanip>
#include <cstdio>

using namespace std;

int main() {
    ios::sync_with_stdio(false);
/*
 redirect stdin and stdout to files for debug ease.
*/
    //freopen("in.txt", "r", stdin);
    //freopen("out.txt", "w", stdout);
    int n, k, i;
    double ci, si, csum, scsum;
    string s;
    bool has, first;
    cin>>n;
    first = true;
    while(n--) {
        cin>>k;
        csum = scsum = 0.0;
        has = false;
        for (i = 0; i < k; ++i) {
            cin>>s>>ci>>si;
            csum += ci;
            scsum += ci * si;
            if (si >= 0 && si < 60) has = true;
        }
        if (first) {
            // nothing
            first = false;
        } else {
            cout<<endl;
        }
        if (has) cout<<"Sorry!"<<endl;
        else {
            cout<<fixed<<setprecision(2)<<scsum / csum<<endl;
        }
    }
    return 0;
}

