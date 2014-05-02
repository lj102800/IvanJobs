#include <cstdio>

double make(int n) {
    double res = 0.0;

    int i;
    for (i = 1; i <= n; ++i)
    {
        if (i % 2 == 1)
            res += 1.0 / i;
        else
            res -=  1.0 / i;
    }
    return res;
}

int main()
{
    int m, n;
    scanf("%d", &m);
    while (m > 0) {
        scanf("%d", &n);
        printf("%.2lf\n",  make(n));
        m--;
    }
    return 0;
}
