// fast long long sqrt, did not get it right now!
long long fastIntSqrt(long long n) {
    long long temp, nHat = 0, b = 0x80000000, bshft = 31;

    do  {
        if (n >= (temp = (((nHat<<1)+b)<<bshft--)))   {
             nHat += b;
             n -= temp;
         }
     } while (b >>= 1);
    return nHat;
}
