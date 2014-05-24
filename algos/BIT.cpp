void build()
{
    for (int i=1;i<=MAX_N;i++)
    {
        BIT[i]=A[i];
        for (int j=i-lowbit(i);j>0;j-=lowbit(j))
            BIT[i]+=BIT[j];
    }
}

void edit(int i, int delta)
{
    for (int j = i; j <= MAX_N; j += lowbit(j))
        BIT[j] += delta;
}

int sum (int k)
{
    int ans = 0;
    for (int i = k; i > 0; i -= lowbit(i))
        ans += BIT[i];
    return ans;
}

