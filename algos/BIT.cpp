/**
 * BIT, Binary Index Tree, so-called Tree Array.
 * we store sequence in A, index starts from 1
 * and we make another array C, that C[i] = sum(C[i - 2^k + 1] to C[i]) 
 * and 2^k = lowbit(i), we make C[0] = 0; for consistent use of sum(a) - sum(b - 1) 
 */
int lowbit(int t)
{
    return t & (-t);
}

void build()
{
    BIT[0] = 0;
    for (int i=1;i<=MAX_N;i++)
    {
        BIT[i]=A[i];
        for (int j=i-1;j>i - lowbit(i) + 1;j--)
            BIT[i]+=A[j];
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

