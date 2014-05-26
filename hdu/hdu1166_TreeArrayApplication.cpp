#include <cstdio>
#include <cmath>
#define MAXN 50100
int n;
int BIT[MAXN];
int A[MAXN];
#define Log(fmt, ...) do {if(false) printf(fmt"\n", ##__VA_ARGS__);} while(false)

int lowbit(int t)
{
    return t & (-t);
}

void build()
{
    BIT[0] = 0;
    for (int i=1;i<=n;i++)
    {
        BIT[i]=A[i];
        for (int j=i-1;j>=i - lowbit(i) + 1;j--) {
        Log("j:%d", j);
            BIT[i]+=A[j];
    }
    Log("BIT[%d]:%d", i, BIT[i]);
    }
}


void edit(int i, int delta)
{
    for (int j = i; j <= n; j += lowbit(j))
        BIT[j] += delta;
}


int sum (int k)
{
    int ans = 0;
    for (int i = k; i > 0; i -= lowbit(i))
        ans += BIT[i];
    return ans;
}

int main() {
    int T, tn;
    scanf("%d", &T);
    for (tn = 1; tn <= T; ++tn) {
        scanf("%d", &n);
        int i, a, b;
        for (i = 1; i <= n; ++i) scanf("%d", &A[i]);
        char cmd[10];
        build();
        printf("Case %d:\n", tn);
        while(scanf("%s", cmd) && cmd[0] != 'E') {
            switch(cmd[0]) {
                case 'A':
                    scanf("%d%d", &a, &b);
                    edit(a, b);        
                break;
                case 'S':
                    scanf("%d%d", &a, &b);
                    edit(a, -b);
                break;
                case 'Q':
                    scanf("%d%d", &a, &b);
                    printf("%d\n", sum(b) - sum(a - 1));
                break;
            }
        }
    }
    return 0;
}
