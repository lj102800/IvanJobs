#include <iostream>
#include <vector>
#include <cstring>
using namespace std;

const int _max = 60;
int c1[_max], c2[_max]; // c1 多项式展开每个阶段的结果， c2 保存c1 和 当前子表达式相乘的结果。
vector<int> cnt, v;
int main()
{
    int t, i, j, k;
    cin>>t;
    while(t--) {
        memset(c1, 0, sizeof(c1)); memset(c2, 0, sizeof(c2));
        cnt.clear(); v.clear();
        for (i = 1; i <= 26; ++i) {
            int tmp;
            cin>>tmp;
            if (tmp != 0) {
                cnt.push_back(tmp); v.push_back(i); // 对个数为0项的过滤。
            }
        }
        if (cnt.empty()) cout<<0<<endl;
        else {
            for (i = 0; i <= 50 && i / v[0] <= cnt[0]; i += v[0]) c1[i] = 1;// 初始化第一个表达式到c1
            int len = (int)cnt.size();
            for (i = 1; i < len; ++i) { // 每一次循环都是将当前c1 和 当前表达式相乘。
                // calc the second expression
                for (j = 0; j <= 50; ++j) { // 遍历c1的所有项
                    for (k = 0; k + j <= 50 && k/v[i] <= cnt[i]; k += v[i]) c2[k+j] += c1[j];// 和当前表达式i相乘。
                }
                for(int j=0; j<=50; j++){// c1 和 c2 交换，保持循环不变式的状态。
                        c1[j] = c2[j];
                        c2[j] = 0;
                }
            }
            int res = 0;
            for (i = 1; i <= 50; ++i) res += c1[i];
            cout<<res<<endl;
        }
    }

    return 0;
}
