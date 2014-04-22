#include <iostream>
#include <string>
#include <stack>
#include <algorithm>
#include <deque>
#include <vector>
#include <cstdio>

using namespace std;
#define L(fmt,...) do {if(false) printf(fmt"\n", ##__VA_ARGS__);}while(false)

bool can(string& src, string& tar) {
    if (src.size() != tar.size()) return false;
    vector<char> va, vb;
    int i;
    int len = src.size();
    for (i = 0; i < len; ++i){
        va.push_back(src[i]);
        vb.push_back(tar[i]);
    }
    sort(va.begin(), va.end());
    sort(vb.begin(), vb.end());
    for (i = 0; i < len; ++i)
        if (va[i] != vb[i]) return false;
    return true;
}
string src, tar, res;
vector<string> vr;
void load(deque<char>& r) {
    int len = r.size();
    int i;
    string re;
    for (i = 0; i < len; ++i){
        re.push_back(r[i]);
        re.push_back(' ');
    }
    L("load:%s", re.c_str());
    vr.push_back(re);
}
void display() {
    sort(vr.begin(), vr.end());
    int len = vr.size();
    int i;
    for (i = 0; i < len; ++i)
        cout<<vr[i]<<endl;
}

void bt(int i, int j, stack<char>& c, deque<char>& q) {
    L("bt(%d, %d)", i, j);
    if (i == (int)src.size() && j == (int)tar.size()) {
        // find one
        load(q);
        return ;
    }

    // do i
    if (i < (int)src.size()) {
        c.push(src[i]);
        q.push_back('i');
        bt(i + 1, j, c, q);
        c.pop();
        q.pop_back();
    }
    // do o
    if (!c.empty()) {
        char v = c.top();

        if (v == tar[j]) {
            q.push_back('o');
            c.pop();
            bt(i, j + 1, c, q);
            c.push(v);
            q.pop_back();
        }
    }
}

int main() {
    ios::sync_with_stdio(false);

    while(cin>>src>>tar) {
        if (!can(src, tar)) {
            cout<<"["<<endl<<"]"<<endl;
            continue;
        }
        stack<char> c;
        deque<char> q;
        vr.clear();
        bt(0, 0, c, q);
        cout<<"["<<endl;
        display();
        cout<<"]"<<endl;
    }
    return 0;
}
