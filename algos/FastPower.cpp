int pow2(int a, int b) {
    int r = 1, base = a;
    while(b != 0) {
        if (b % 2) r *= base;
        base *= base;
        b /= 2;
    }
    return r;
}

int pow3(int a, int b) {
    int r = 1, base = a;
    while(b != 0) {
        if (b & 1) r *= base;
        base *= base;
        b >>= 1;
    }
    return r;
}

int pow4(int x, int n) {
    int result;
    if (n == 0) return 1;
    else {
        while((n & 1) == 0) {
            n >>= 1;
            x *= x;
        }
    }
    result = x;
    n >>= 1;
    while(n != 0) {
        x *= x;
        if ((n & 1) != 0) result *= x;
        n >>= 1;
    }
    return result;
}
