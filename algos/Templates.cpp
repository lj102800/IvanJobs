第一章 常用函数和STL
一， 常用函数
#include <cstdio>
int getchar(); // 读取一个字符，一般用来去掉无用字符。
char* gets(char* str); // 读取一行字符串

#include <cstdlib>
void* malloc(size_t size); // 动态内存分配， 开辟大小为size的空间
void qsort(void* buf, size_t num, size_t size, int(*compare)(const void*, const void*)); // 快速排序
Sample.
int compare_ints(const void* a, const void* b) {
	int* arg1 (int*)a; int* arg2 = (int*)b;
	if (*arg1 < *arg2) return -1;
	else if (*arg1 == *arg2) return 0;
	else return 1;
}
int array[] = {-2, 99, 0, -743, 2, 3, 4}; int array_size = 7;
qsort(array, array_size, sizeof(int), compare_ints);

#include <cmath>
// 求反正弦， arg [-1, 1], 返回值 [-pi/2, pi/2]
double asin(double arg);
double sin(double arg);
double exp(double arg);
double log(double num);
double sqrt(double num);
double pow(double base, double exp);

#include <cstring>
void* memset(void* buffer, int ch, size_t count);
memset(the_array, 0, sizeof(the_array));
int sprint(char* buffer, const char* format, ...);
sprintf(s, "%d%d", 123, 4567);
int sscanf(const char* buffer, const char* format, ...);
Sample:
char result[100] = "24 hello", str[100]; int num;
sprintf(result, "%d %s", num , str);
int strcmp(const char* str1, const char* str2);

二，常用STL
[标准container概要]
vector<T>	大小可变的向量， 类似数组里的用法。
list<T>		双向链表
queue<T>	队列， empty(), front(), pop(), push()
stack<T>	栈， empty(), top(), pop(), push()
priority_queue<T>	优先队列
set<T>		集合
map<key, val>	关联数组，常用来作为hash映射

[标准algorithm摘录]
for_each()	对每一个元素都调用一个函数
find()
replace()	用新的值替换元素， O(N)
copy()		复制元素， O(N)
remove()	移除元素
reverse()	倒置元素
sort()		排序， O(NlogN)
partial_sort()	部分排序
binary_search() 二分查找
merge()		合并有序的序列，O(N)

[C++ String 摘录]
copy()	从别的字符串拷贝
empty()	判断字符串是否为空
erase() 从字符串移除元素
find() 查找元素
insert() 插入元素
length() 字符串长度
replace() 替换元素
substr() 取子字符串
swap()	交换字符串

第二章， 重要公式和定理
1. Fibonacci Number
0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610...
Formula:
F(0) = 0
F(1) = 1
F(i) = F(i - 1) + F(i - 2)
F(n) = ((1 + sqrt(5))^n - (1 - sqrt(5))^n)/(2^n * sqrt(5)) = [1/sqrt(5) * (((1 + sqrt(5))/2)^n)]
2. Lucas Number
1, 3, 4, 7, 11, 18, 29, 47, 76, 123 ...
Formula:
L(n) = ((1 + sqrt(5))/2)^n + ((1 - sqrt(5))/2)^n
3. Catalan Number
1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796, 58786, 208012 ...
Formula:
Cat(n) = C(2n, n)/(n + 1)
Cat(n) = sum{Cat(i) * Cat(n - 1 - i) | i = 0 -> n-1}
Application:
1) 将n + 2边形沿弦切割成n 个三角形的不同切割数。
2）n + 1 个数相乘， 给每两个元素加上括号的不同方法数
3) n 个节点的不同形状的二叉树数
4) 从n * n 方格的左上角移动到右下角不升路径数。
4. Stirling Number(Second Kind)
S(n, m)表示含有n个元素的集合划分成m个集合的情况数
或者是n个有标号的球，放到m个无标号的盒子中，要求无一为空，其不同的方案数
Formula:
S(n, m):
if m = 0 || n < m, then S(n, m) = 0
if n > m >= 1 , then S(n, m) = S(n - 1, m - 1) + m * S(n - 1, m)

S(n, m) = 1/(m!) * sum{(-1)^i * C(m, i) * (m - i)^n}

Special Cases:
S(n, 0) = 0
S(n, 1) = 1
S(n, 2) = 2^(n - 1) - 1
S(n, 3) = (1/6) * (3^n - 3 * 2^n + 3)
S(n, n - 1) = C(n, 2)
S(n, n) = 1

5. Bell Number
n 个元素集合所有的划分数
Formula:
B(n) = sum{S(n, i) | i = 0 -> n}

6. Stirling's Approximation
n! = sqrt(2*pi*n) * (n/e)^n

7. Sum of Reciprocal Approximation
EulerGamma = 0.57721566490153286060651209
sum{1/i | i = 1 -> n} = ln(n) + EulerGamma(n->inf)

8. Young Tableau
Young Tableau(杨氏图表) 是一个矩阵， 它满足条件：
如果格子[i, j]没有元素， 则[i + 1, j]也一定没有元素
如果格子[i, j]有元素a[i, j], 则[i + 1, j]要么没有元素，要么a[i + 1, j] > a[i, j]
Y[n] 代表n个数所组成的杨氏图表的个数。

Formula:
Y(1) = 1
Y(2) = 2
Y(n) = Y(n - 1) + (n - 1) * Y(n - 2) (n > 2)

9. 整数划分
将整数n分成k份， 且每份不能为空，任意两种分法不能相同
1) 不考虑顺序
for (int p = 1; p <=n; p++)
	for (int i = p; i <= n; i++)
		for (int j = k; j >= 1; j--)
			dp[i][j] += dp[i - p][j - 1];
cout<<dp[n][k]<<endl;
2) 考虑顺序
dp[i][j] = dp[i - k][ j - 1]; (k = 1...i)
3) 若分解出来的每个数均有一个上限m
dp[i][j] = dp[i - k][j - 1]; (k = 1...m)

10. 错排公式
D(1) = 0
D(2) = 1
D(n) = (n - 1) * (D(n - 1) + D(n - 2))

11. 三角形内切圆半径公式
p = (a + b + c) / 2;
s = sqrt(p * (p - a) * (p - b) * (p - c));
r = 2*s / (a + b + c);

12. 三角形外接圆半径公式
R = (abc)/(4s)

13. 圆内接四边形面积公式
p = (a + b + c + d) / 2
s = sqrt((p - a) * (p - b) * (p - c) * (p - d));

14. 基础数论公式
1） 模取幂
a^n % b = (((a % b) * a)%b...)% b
2) n的约数的个数
若 n满足 n = p1^n1 + p2^n2 + ... + pm^nm, 则n的约数的个数为
(n1 + 1)(n2 + 1)...(nm + 1)


第三章 大数模板

typedef int hugeint
// 应不大于， 以防乘法时溢出
const int Base = 1000;
const int Capacity = 1000;

struct xnum {
	int Len;
	int Data[Capacity];
	xnum() : Len(0) {}
	xnum(const xnum& V) : Len(V.Len) {
		memcpy(Data, V.Data, Len * sizeof(*Data));
	}
	xnum(int V) : Len(0) {
		for (; V > 0; V /= Base) Data[Len++] = V % Base;
	}

	xnum(char S[]);
	xnum& operator=( const xnum& V) {
		Len = V.Len;
		memcpy(Data, V.Data, Len * sizeof(*Data));
		return *this;
	}
	int& operator[] (int Index) { return Data[Index];}
	int operator[] (int Index) { return Data[Index];}
	void print() {
		printf("%d", Len == 0 ? 0 : Data[Len - 1]);
		for (int i = Len - 2; i >= 0; i--)
			for (int j = Base/10; j > 0; j /= 10)
				printf("%d", Data[i]/j%10);
	}
};

xnum::xnum(char S[]) {
	int I, J;
	Data[Len = 0] = 0;
	J = 1;
	for (I = strlen(S) - 1; I >= 0; I--) {
		Data[Len] += (S[i] - '0') * J;
		J *= 10;
		if (J >= Base) J = 1, Data[++Len] = 0;
	}
	if (Data[Len] > 0) Len++;
}

int compare(const xnum& A, const xnum& B) {
	int I;
	if (A.Len != B.Len) return A.Len > B.Len ? 1 : -1;
	for (I = A.Len - 1; I >= 0 && A[I] == B[I]; I--);
	if (I < 0) return 0;
	return A[I] > B[I] ? 1 : -1;
}

xnum operator+(const xnum& A, const xnum& B) {
	xnum R;
	int I;
	int Carry = 0;
	for (I = 0; I < A.Len || I < B.Len || Carry > 0; I++) {
		if (I < A.Len) Carry += A[I];
		if (I < B.Len) Carry += B[I];
		R[I] = Carry % Base;
		Carry /= Base;
	}
	R.Len = I;
	return R;
}

xnum operator-(const xnum& A, const xnum& B) {
	xnum R;
	int Carry = 0;
	R.Len = A.Len;
	int I;
	for (I = 0; I < R.Len; I++) {
		R[i] = A[I] - Carry;
		if (I < B.Len) R[I] -= B[I];
		if (R[I] < 0) Carry = 1, R[I] += Base;
		else Carry = 0;
	}
	while(R.Len > 0 && R[R.Len - 1] == 0) R.Len--;
	return R;
}

xnum operator*(const xnum& A, const int B) {
	int I;
	if (B == 0) return 0;
	xnum R;
	hugeint Carry = 0;
	for (I = 0; I < A.Len || Carry > 0; I++) {
		if (I < A.Len) Carry += hugeint(A[I]) * B;
		R[I] = Carry % Base;
		Carry /= Base
	}
	R.Len = I;
	return R;
}
xnum operator*(const xnum& A, const xnum& B) {
	int I;
	if (B.Len == 0) return 0;
	xnum R;
	for (I = 0; I < A.Len; I++) {
		hugeint Carry = 0;
		for (int J = 0; J < B.Len || Carry > 0; J++) {
			if (J < B.Len) Carry += hugeint(A[I]) * B[J];
			if (I + J < R.Len) Carry += R[I + J];
			if (I + J >= R.Len) R[R.Len++] = Carry % Base;
			else R[I + J] = Carry % Base;
			Carry /= Base;
		}
	}
	return R;
}


xnum operator/(const xnum& A, const int B) {
	xnum R;
	int I;
	hugeint C = 0;
	for (I = A.Len - 1; I >= 0; I--){
		C = C * Base + A[I];
		R[I] = C / B;
		C %= B;
	}
	R.Len = A.Len;
	while(R.Len > 0 && R[R.Len - 1] == 0) R.Len--;
	return R;
}

xnum operator/(const xnum& A, const xnum& B) {
	int I;
	xnum R, Carry = 0;
	int Left, Right, Mid;
	for (I = A.Len - 1; I >= 0; I--) {
		Carry = Carry * Base + A[I];
		Left = 0;
		Right = Base - 1;
		while(Left < Right) {
			Mid = (Left + Right + 1) / 2;
			if (compare(B*Mid, Carry) <= 0) Left = Mid;
			else Right = Mid - 1；
		}
		R[I] = Left;
		Carry = Carry - B * Left;
	}
	R.Len = A.Len;
	while(R.Len > 0 && R[R.Len - 1] == 0) R.Len--;
	return R;
}

xnum operator%(const xnum& A, const xnum& B) {
	int I;
	xnum R, Carry = 0;
	int Left, Right, Mid;
	for (I = A.Len - 1; I >= 0; I--) {
		Carry = Carry * Base + A[I];
		Left = 0;
		Right = Base - 1;
		while(Left < Right) {
			Mid = (Left + Right + 1)/2;
			if (compare(B*Mid, Carry) <= 0) Left = Mid;
			else Right = Mid - 1;
		}
		R[I] = Left;
		Carry = Carry - B * Left;
	}
	R.Len = A.Len;
	while(R.Len > 0 && R[R.Len - 1] == 0) R.Len--;
	return Carry;
}

istream& operator>>(istream& In, xnum& V) {
	char Ch;
	for (V = 0; In >>Ch;) {
		V = V * 10 + (Ch - '0');
		if (cin.peek() <= '') break;
	}
	return In;
}

ostream& operator<<(ostream& Out, const xnum& V) {
	int I;
	Out<<(V.Len == 0 ? 0 : V[V.Len - 1]);
	for (I = V.Len - 2; I >= 0; I--)
		for (J = Base / 10; J > 0; J /= 10)
			Out<<V[I]/J%10;
	return Out;
}

xnum gcd(xnum a, xnum b) {
	if (compare(b, 0) == 0) return a;
	else return gcd(b, a % b);
}

int div(char* A, int B) {
	int I;
	int C = 0;
	int Alen = strlen(A);
	for (I = 0; I < Alen; I++) {
		C = C * Base + A[I] - '0';
		C %= B;
	}
	return C;
}

xnum C(int n, int m) {
	int i;
	xnum sum = 1;
	for (i = n; i >= n - m + 1; i--)
		sum = sum * i;
	for (i = 1; i <= m; i++)
		sum = sum / i;
	return sum;
}



#define MAXN 9999
#define DLEN 4
class BigNum {
	private:
		int a[1000];
		int len;
	public:
		BigNum() {len = 1; memset(a, 0, sizeof(a));}
		BigNum(const int);
		BigNum(const char*);
		BigNum(const BigNum&);
		BigNum& operator=(const BigNum&);
		BigNum operator+(const BigNum&) const;
		BigNum operator-(const BigNum&) const;
		BigNum operator*(const BigNum&) const;
		BigNum operator/(const BigNum&) const;
		BigNum operator^(const int&) const;
		int operator%(const int &) const;
		bool operator>(const BigNum& T) const;
		void print();
};

BigNum::BigNum(const int b) {
	int c, d = b;
	len = 0;
	memset(a, 0, sizeof(a));
	while(d > MAXN) {
		c = d - ((d / (MAXN + 1)) * (MAXN + 1));
		d = d / (MAXN + 1); a[len++] = c;
	}
	a[len++] = d;
}

BigNum::BigNum(const char* s) {
	int t, k, index, I, i;
	memset(a, 0, sizeof(a));
	I = strlen(s);
	len = I/DLEN;
	if (I % DLEN) len++;
	index = 0;
	for (i = I - 1; i >= 0; i -= DLEN) {
		t = 0; k = i - DLEN + 1;
		if (k < 0) k = 0;
		for (int j = k; j <= i; j++)
			t = t * 10 + s[j] - '0';
		a[index++] = t;
	}
}

BigNum::BigNum(const BigNum& T) : len(T.len) {
	int i;
	memset(a, 0, sizeof(a));
	for(i = 0; i < len; ++i) a[i] = T.a[i];
}

BigNum& BigNum::operator=(const BigNum& n) {
	len = n.len;
	memset(a, 0, sizeof(a));
	int i;
	for (i = 0; i < len; ++i)
		a[i] = n.a[i];
	return *this;
}

BigNum BigNum::operator+(const BigNum& T) const {
	BigNum t(*this);
	int i, big;
	big = T.len > len ? T.len : len;
	for (i = 0; i < big; i++) {
		t.a[i] += T.a[i];
		if (t.a[i] > MAXN) {
			t.a[i + 1]++;
			t.a[i] -= MAXN + 1;
		}
	}
	if (t.a[big] != 0) t.len = big + 1;
	else t.len = big;
	return t;
}

BigNum BigNum::operator-(const BigNum& T) const {
	int i, j, big;
	bool flag;
	BigNum t1, t2;
	if (*this > T) {
		t1 = *this;
		t2 = T;
		flag = 0;
	} else {
		t1 = T;
		t2 = *this;
		flag = 1;
	}
	big = t1.len;
	for (i = 0; i < big; ++i) {
		if (t1.a[i] < t2.a[i]) {
			j = i + 1;
			while(t1.a[j] == 0) j++;
			t1.a[j--]--;
			while(j > i) t1.a[j--] += MAXN;
			t1.a[i] += MAXN + 1 - t2.a[i];
		} else t1.a[i] -= t2.a[i];
	}
	t1.len = big;
	while(t1.a[len - 1] == 0 && t1.len > 1) {
		t1.len--;
		big--;
	}
	if (flag) t1.a[big - 1] = 0 - t1.a[big - 1];
	return t1;
}

BigNum BigNum::operator*(const BigNum& T) const {
	BigNum ret;
	int i, j, up;
	int temp, temp1;
	for (i = 0; i < len; ++i) {
		up = 0;
		for (j = 0; j < T.len; ++j) {
			temp = a[i] * T.a[j] + ret.a[i + j] + up;
			if(temp > MAXN) {
				temp1 = temp - temp /(MAXN + 1) * (MAXN + 1);
				up = temp / (MAXN + 1);
				ret.a[i + j] = temp1;
			} else {
				up = 0;
				ret.a[i + j] = temp;
			}
		}
		if (up != 0)
			ret.a[i + j] = up;
	}
	ret.len = i + j;
	while(ret.a[ret.len - 1] == 0 && ret.len > 1) ret.len--;
	return ret;
}

BigNum BigNum::operator/(const int& b) const {
	BigNum ret;
	int i, down = 0;
	for (i = len - 1; i >= 0; i--) {
		ret.a[i] = (a[i] + down * (MAXN + 1)) / b;
		down = a[i] + down * (MAXN + 1) - ret.a[i] * b;
	}
	ret.len = len;
	while(ret.a[ret.len - 1] == 0 && ret.len > 1) ret.len--;
	return ret;
}

int BigNum::operator %(const int & b) const {
	int i, d = 0;
	for (i = len - 1; i >= 0; i--) {
		d = ((d * (MAXN + 1)) % b + a[i]) % b;
	}
	return d;
}

BigNum BigNum::operator^(const int & n) const {
	BigNum t, ret(1);
	if (n < 0) exit(-1);
	if (n == 0) return 1;
	if (n == 1) return *this;
	int m = n;
	while(m > 1) {
		t = *this;
		int i;
		for (i = 1; i<<1 <= m; i <<= 1) {
			t = t*t;
		}
		m -= i;
		ret = ret * t;
		if (m == 1) ret = ret * (*this);
	}
	return ret;
}

bool BigNum::operator>(const BigNum& T) const {
	int ln;
	if (len > T.len) return true;
	else if (len == T.len) {
		ln = len - 1;
		while(a[ln] == T.a[ln] && ln >= 0) ln--;
		if (ln >= 0 && a[ln] > T.a[ln]) return true;
		else return false;
	} else return false;
}

void BigNum::print() {

	int i;
	cout<<a[len - 1];
	for (i = len - 2; i >= 0; i--) {
		cout.width(DLEN);
		cout.fill('0');
		cout<<a[i];
	}
}

// 读取整数
const int ok = 1;
int get_val(int& ret) {
	ret = 0;
	char ch;
	while((ch = getchar()) > '9' || ch < '0');
	do {
		ret = ret * 10 + ch - '0';
	} while((ch = getchar()) <= '9' && ch >= '0');
	return ok;
}

int get_val(int & ret) {
	ret = 0;
	char ch;
	bool neg = false;
	while(((ch=getchar()) > '9' || ch < '0') && ch != '-');
	if (ch == '-') {
		neg = true;
		while((ch = getchar()) > '9' || ch < '0') ;
	}
	do {
		ret = ret * 10 + ch - '0';
	} while((ch = getchar()) <= '9' && ch >= '0');
	ret = (neg ? -ret : ret);
	return ok;
}

// 读取整数，可判断EOF和EOL
const int eof = -1;
const int eol = -2;

int get_val(int & ret) {
	ret = 0;
	char ch;
	while(((ch = getchar()) > '9' || ch < '0') && ch != EOF);
	if (ch == EOF) return eof;
	do {
		ret = ret * 10 + ch - '0';
	} while((ch = getchar()) <= '9' && ch >= '0');
	if (ch == '\n') return eol;
	return ok;
}

// 读取浮点数
int get_val(double & ret) {
	ret = 0;
	double base = 0.1;
	char ch;
	bool dot = false, neg = false;
	while(((ch = getchar()) > '9' || ch < '0') && ch != '.' && ch != '-');
	if (ch == '-') {
		neg = true;
		while(((ch = getchar()) > '9' || ch < '0') && ch != '.' && ch != '-');
	}
	do {
		if (ch == '.' ){
			dot = true;
			continue;
		}
		if (dot) {
			ret += (ch - '0') * base;
			base *= 0.1;
		} else ret = ret * 10 + (ch - '0');
	} while(((ch = getchar()) <= '9' && ch >= '0') || ch == '.');
	ret = (neg ? -ret : ret);
	return ok;
}

typedef long long LL;
LL MultiMod(LL a, LL b, LL c) {
	LL ret = 0, d = a;
	for (; b; b >>= 1, d <<= 1, d %= c) if (b & 1) ret = (ret + d) % c;
	return ret;
}

// 128-bits integer's power with mod in O(64 * LogN)
LL ModPower(LL base, LL exp, LL mod) {
	LL ret = 1;
	for (; exp; exp >>= 1, base = MultiMod(base, base, mod)) if (exp & 1) ret = MultiMod(ret, base, mod);
	return ret;
}


第四章，数论算法
1. Greatest Common Divisor 最大公约数
int GCD(int x, int y) {
	int t;
	while(y > 0) {
		t = x % y;
		x = y;
		y = t;
	}
	return x;
}

2. Prime 素数判断
bool is_prime(int u) {
	if (u == 0 || u == 1) return false;
	if (u == 2) return true;
	if (u % 2 == 0) return false;
	for (int i = 3; i <= sqrt(u); i += 2) if (u % i == 0) return false;
	return true;
}
3. Sieve Prime 素数筛法
const int M = 1000; // M : size
bool mark[M]; // true : prime number
void sieve_prime() {
	memset(mark, true, sizeof(mark));
	mark[0] = mark[1] = false;
	for (int i = 2; i <= sqrt(M); i++) if (mark[i]) {
		for (int j = i * i; j < M; j += i) mark[j] = false;
	}
}

4. Module Inverse 模逆元
int Inv(int a, int n) {
	int d, x, y;
	d = extended_euclid(a, n, x, y);
	if (d == 1) return (x % n +n) % n;
	else return -1; // no solution
}

5. Extened Euclid 扩展欧几里得算法
// 如果GCD(a, b) = d, 则存在 x, y, 使得 d = ax + by
// extended_euclid(a, b) = ax + by
int extended_euclid(int a, int b, int& x, int &y) {
	int d;
	if (b == 0) { x = 1; y = 0; return a;}
	d = extended_euclid(b, a % b, y, x);
	y -= a / b * x;
	return d;
}

6. Modular Linear Equation 模线性方程（同余方程）
// 如果GCD(a, b) 不能整除c, 则ax + by = c 没有整数解
void modular_linear_equation(int a, int b, int n) {
	int d, x, y, x0, gcd;
	gcd = GCD(a, n);
	if (b % gcd != 0) {
		cout<<"no solution"<<endl;
		return ;
	}
	a /= gcd; b /= gcd; n /= gcd;
	d = extened_euclid(a, n, x, y);
	if (b % d == 0) {
		x0 = (x * (b / d)) % n; // x0 : basic solution
		int ans = n; // min x = (x0 % (n / d) + (n / d) % (n / d))
		for (int i = 0; i < d; i++) {
			ans = (x0 + i * (n / d)) % n;
			cout<<ans<<endl;
		}
	} else cout<<"no solution"<<endl;
}

7. Chinese Remainder Theorem 中国余数定理
int chinese_remainder(int b[], int w[], int len) {
	int i, d, x, y, m, n, ret;
	ret = 0; n = 1;
	for (i = 0; i < len; ++i) n *= w[i];
	for (i = 0; i < len; ++i) {
		m = n / w[i];
		d = extended_euclid(w[i], m, x, y);
		ret = (ret + y * m * b[i]) % n;
	}
	return (n + ret%n) % n;
}

// Pku 2891 Strange Way to Express Integers
LL chinese_remainder2() {
	int i, j;
	if (n == 1) return r[0];
	LL m, x, apre;
	x = modular_linear_equation(a[0], r[1] - r[0], a[1]);
	if(x == -1) return -1;
	m = x * a[0] + r[0];
	apre = LCM(a[0], a[1]);
	for(i = 2; i < n; i++) {
		x = modular_linear_equation(apre, r[i] - m, a[i]);
		if (x == -1) return -1;
		m = x * apre + m;
		apre = LCM(apre, a[i]);
	}
	return m;
}

8. Euler Function 欧拉函数
// 求 1...n-1 中与n互质的个数
int euler(int n) {
	int ans = 1;
	int i;
	for (i = 2; i*i <= n; i++) {
		if (n % i == 0) {
			n /= i;
			ans *= i - 1;
			while(n % i == 0) {
				n /= i;
				ans *= i;
			}
		}
	}
	if (n > 1) ans *= n - 1;
	return ans;
}

9. Farey 总数
// 求MAX以内所有Farey 的总数
const int MAX = 1000100;
int n;
bool num[1100]; // sqrt(MAX)
int prime[1100], total;
__int64 f[MAX], inc[MAX];

void cal_prime() {
	int i, j;
	memset(num, false, sizeof(num));
	total = 0;
	for (i = 2; i < 1100; i++) {
		if (!num[i]) {
			prime[total++] = i;
			j = i + i;
			while(j < 1100) {
				num[j] = true;
				j += i;
			}
		}
	}
}

void cal_farey() {
	int i, j, k;
	inc[1] = 1;
	for (i = 2; i < MAX; i++) {
		for (j = 0; j < total; j++) {
			if(i % prime[j] == 0) {
				k = i / prime[j];
				if (k % prime[j] == 0) inc[i] = inc[k] * prime[j];
				else inc[i] = inc[k] * (prime[j] - 1);
				break;
			}
		}
		if (j == total) inc[i] = i - 1;
	}
	f[1] = 0;
	for (i = 2; i < MAX; i++) f[i] = f[i - 1] + inc[i];
}

int main() {
	cal_prime();
	cal_farey();
	while(scanf("%d", &n) , n) {
		printf("%I64d\n", f[n]);
	}
	return 0;
}


10. Farey 序列构造
// 构造5000 以内的Farey序列
const int MAX = 8000000;
int total;
int n, k;
int farey[2][MAX];
void make_farey_seq(int x1, int y1, int x2, int y2) {
	if (x1 + x2 > n || y1 + y2 > n) return;
	make_farey_seq(x1, y1, x1 + x2, y1 + y2);
	total++;
	farey[0][total] = x1 + x2;
	farey[1][total] = y1 + y2;
	make_farey_seq(x1 + x2, y1 + y2, x2, y2);
}
int main() {
	int t;
	scanf("%d%d", &n, &t);
	total = 1;
	farey[0][1] = 0;
	farey[1][1] = 1;
	make_farey_seq(0, 1, 1, 1);
	farey[0][total + 1] = 1;
	farey[1][total + 1] = 1;
	total++;
	while(t--) {
		scanf("%d", &k);
		if (k > total) puts("No Solution");
		else printf("%d/%d\n", farey[0][k], farey[1][k]);
	}
}

11. Miller_Rabbin 素数测试， Pollard_rho 因式分解
typedef __int64 I64
const char* pformat = "%I64d";
I64 big_rand(I64 m) {
	I64 x = rand();
	x *= rand();
	if (x < 0) x = -x;
	return x %= m;
}

// x *y % n
I64 mod_mul(I64 x, I64 y, I64 n) {
	if (x == 0 || y == 0) return 0;
	return (((x&1)*y)%n +(mod_mul(x>>1, y, n)<<1)%n)%n;
}

// x^y % n
I64 mod_exp(I64 x, I64 y, I64 n) {
	I64 ret = 1;
	while(y) {
		if (y & 1) ret = mod_mul(ret, x, n);
		x = mod_mul(x, x, n);
		y >>= 1;
	}
	return ret;
}

bool Miller_Rabbin(I64 n) { // O(times * (logN)^3)
	I64 i, j, x, m, k;
	if (n == 2) return true;
	if (n < 2 || !(n&1)) return false;
	m = n - 1; k = 0;
	while(!(m&1)) m>>=1, k++; // binary scan
	for (i = 0; i < 4; ++i) { // test times
		x = big_rand(n - 2) + 2;
		x = mod_exp(x, m, n);
		if (x == 1) continue;
		for (j = 0; j < k; j++) {
			if (x == n - 1) break;
			x = mod_mul(x, x, n);
		}
		if (j >= k) return false;
	}
	return true;
}


I64 gcd(I64 x, I64 y) {
	if (x > y) std::swap(x, y);
	while(x) {
		I64 t = y % x;
		y = x;
		x = t;
	}
	return y;
}

I64 func(I64 x, I64 m) {
	return (mod_mul(x, x, m) + 1) % m;
}

I64 Pollard(I64 n) {
	if (Miller_Rabbin(n)) return n;
	if (!(n&1)) return 2;
	I64 i, x, y, ret;
	i = 1;
	while(true) {
		x = i++;
		y = func(x, n);
		ret = gcd(y - x, n);
		while(ret == 1) {
			x = func(x, n);
			y = func(func(y, n), n);
			ret = gcd((y - x + n) % n, n) % n;
		}
		if (0 < ret && ret < n) return ret;
	}
}

I64 factor[100], nfac, minfac;

void cal_factor(I64 n) {
	I64 x = Pollard(n);
	if (x == n) { minfac = min(minfac, x); return ;}
	cal_factor(x);
	cal_factor(n / x);
}

void print_factor(I64 n) {
	I64 i;
	nfac = 0;
	cal_factor(n);
	std::sort(factor, factor + nfac);
	for (i = 0; i < nfac; i++) {
		if (i > 0) putchar(' ');
		printf(pformat, factor[i]);
	}
	puts("");
}

const I64 lim = 100000;

int main() {
	I64 n, t, i;
	srand((unsigned)time(NULL));
	scanf(pformat, &t);
	while(t--) {
		scanf(pformat, &n);
		if (Miller_Rabbin(n)) puts("Prime");
		else {
			if (!(n&1)) puts("2");
			else {
				for (minfac = 3; minfac < lim && n % minfac ; minfac += 2);
				if (minfac >= lim) {
					I64 rm = sqrt(1.0 * n);
					if (rn * rn == n) {
						minfac = rn;
						cal_factor(rn);
					} else {
						minfac = n;
						cal_factor(n);
					}
				}
				printf(pformat, minfac);
				puts("");
			}
		}
	}
}

第五章 图论算法
1. 最小生成树 （Kruscal算法）
#include <iostream>
#include <algorithm>
#include <cstdio>
#include <cmath>
using namespace std;
struct struct_edges {
	int bv, tv;
	double w;
};

struct_edges edges[10100];//边集
struct struct_a {
	double x, y;
};
struct_a arr_xy[101];
int point[101], n, e; // n 顶点数， e 边数（注意是无向网络）
double sum;
int kruscal_f1(int point[], int v) {
	int i = v;
	while(point[i] > 0) i = point[i];
	return i;
}

bool UDlesser(struct_edges a, struct_edges b) {
	return a.w < b.w;
}

void kruscal(){ // 只需要准备好n, e, 递增的边集edges[]即可使用
	int v1, v2, i, j;
	for (i = 0; i < n; i++) point[i] = 0;
	i = j = 0;
	while(j < n - 1 && i < e) {
		v1 = kruscal_f1(point, edges[i].bv);
		v2 = kruscal_f1(point, edges[i].tv);
		if (v1 != v2) {
			sum += edges[i].w;
			point[v1] = v2;
			j++;
		}
		i++;
	}
}

int main() {
	int k, i, j;
	cin>>n;
	k = 0;
	while(n != 0) {
		sum = 0;
		k++;
		for (i = 0; i < n; i++)
            cin>>arr_xy[i].x>>arr_xy[i].y;
        e = 0;
        for (i = 0; i < n; i++)
            for (j = i + 1; j < n; j++) {
                if (i == j) continue;
                edges[e].bv = i;
                edges[e].tv = j;
                edges[e].w = sqrt((arr_xy[i].x - arr_xy[j].x) * (arr_xy[i].x - arr_xy[j].x) +
                                  (arr_xy[i].y - arr_xy[j].y) * (arr_xy[i].y - arr_xy[j].y));
                e++;
            }
            sort(edges, edges + e, UDlesser);
            kruscal();
            printf("Case #%d:\n", k);
            printf("The minimal distance is :%.2f\n", sum);
            cin>>n;
            if (n != 0) printf("\n");
	}
	return 0;
}

2. 最小生成树 （Prim算法）
#include <iostream>
#include <cmath>
#include <cstdio>
using namespace std;

double sum, arr_list[101][101], min;
int i, j, k = 0, n;
struct struct_a {
    float x, y;
};
struct_a arr_xy[101];
struct struct_b {
    int point;
    float lowcost;
};
struct_b closedge[101];

void prim(int n) {//prim需要准备：n顶点数， arr_list[][]顶点的邻接矩阵也是从0开始计数的
    int i, j, k;
    k = 0;
    for (j = 0; j < n; j++) {
        if (j != k) {
            closedge[j].point = k;
            closedge[j].lowcost = arr_list[k][j];
        }
    }
    closedge[k].lowcost = 0;
    for (i = 0; i < n; ++i) {
        min = 10000;
        for (j = 0; j < n; j++) {
            if (closedge[j].lowcost != 0 && closedge[j].lowcost < min) {
                k = j;
                min = closedge[j].lowcost;
            }
        }
        sum += closedge[k].lowcost; // 不要改成 sum += min; sum 即为所求值
        closedge[k].lowcost = 0;
        for (j = 0; j < n; j++) {
            if (arr_list[k][j] < closedge[j].lowcost) {
                closedge[j].point = k;
                closedge[j].lowcost = arr_list[k][j];
            }
        }
    }
}

int main() {
    cin>>n;
    while(n != 0) {
        sum = 0;
        k++;
        for (i = 0; i < n; ++i)
            cin>>arr_xy[i].x>>arr_xy[i].y;
        for (i = 0; i <n; i++)
            for (j = 0; j < n; j++)
                arr_list[i][j] = arr_list[j][i] = sqrt((arr_xy[i].x - arr_xy[j].x)*(arr_xy[i].x - arr_xy[j].x) +
                                                       (arr_xy[i].y - arr_xy[j].y)*(arr_xy[i].y - arr_xy[j].y));
        prim(n);
        cout<<"Case #"<<k<<":"<<endl;
        printf("The minimal distance is:%.2f\n", sum);
        cin>>n;
        if (n != 0) printf("\n");
    }
    return 0;
}

3. 单源最短路径(Bellman-ford算法)
struct node {
    int e, v;
    node(int a = 0, int b = 0):e(a), v(b) {}
};
vector< vector<node> >path;
int n , p, q;
int dist[1000100];
/*******************************
* SPFA (Shortest Path Faster Algorithm)
* Bellman-Ford 算法的一种队列实现， 减少了不必要的冗余计算
* 返回值false， 说明队列不为空， 存在负权环
*********************************/
bool SPFA() {
    int i, j, k, now, l;
    node next;
    bitset<1000100> vis;
    queue<int> SQ;
    memset(dist, -1, sizeof(dist));
    SQ.push(1);
    vis[1] = true;
    dist[1] = 0;

    for (i = 0; i <= n; i++) {
        l = SQ.size();
        if (l == 0) break;
        while(l--) {
            now = SQ.front();
            SQ.pop();
            vis[now] = false;
            for (j = path[now].size() - 1; j >= 0; j--) {
                next = path[now][j];
                if (dist[next.e] == -1 || dist[next.e] > dist[now] + next.v) {
                    dist[next.e] = dist[now] + next.v;
                    if (!vis[next.e]) {
                        SQ.push(next.e);
                        vis[next.e] = true;
                    }
                }
            }
        }
    }
    return SQ.empty();
}

4. 单源最短路径 （Dijkstra算法）
/************************
 贪心， 不能有负权。
***************************/
int matrix[200][200], n;

void Dijkstra(int x, int y) {// 起点 Vx, 终点 Vy
    int i, j, k, path[40000], mark[40000];
    int min, dist[40000];
    for (i = 1; i <= n; i++) {
        mark[i] = 0;
        dist[i] = matrix[x][i];
        path[i] = x;
    }
    mark[x] = 1;
    do {
        min = 30000;
        k = 0;
        for (i = 1; i <= n; i++) if (mark[i] == 0 && dist[i] < min) {
            min = dist[i];
            k = i;
        }
        if (k) {
            mark[k] = 1;
            for (i = 1; i <= n; i++) if (matrix[k][i] < 30000 && min +matrix[k][i] < dist[i]) {
                dist[i] = min + matrix[k][i];
                path[i] = k;
            }
        }
    } while(k);
    cout<<dist[y]<<endl;
    // 如果希望得到路径，加入如下代码：
    do {
        cout<<k<<"<--";
        k = path[k];
    } while(k != x);
    cout<<x<<endl;
}

5. 全源最短路径(Folyd算法)
/*****************************
DP, O(N^3)
******************************/
// 初始化
// path[i][j] = j;
void Floyd() {
    int i, j, k;
    for(k = 0; k < vertex_number; k++) {
        for (i = 0; i < vertex_number; i++) {
            for (j = 0; j < vertex_number;j++) {
                if ((graph[i][k] == -1) || (graph[k][j] == -1)) continue;
                if ((graph[i][j] == -1) || (graph[i][j] > graph[i][k] + graph[k][j])) {
                    graph[i][j] = graph[i][k] + graph[k][j]; /* 最短路径值 */
                    path[i][j] = k; /* 最短路径 */
                }
            }
        }
    }
}


6. 拓扑排序
// degree[] 每个节点的入度
// f[]  每个节点所在的层

void Toplogical_sort() {
    int i, j;
    bool p = true;
    top = 0;
    while(p) {
        p = false;
        top++;
        for (i = 1; i <= n; i++) if (degree[i] == 0) {
            p = true;
            f[i] = top;
        }
        for (i = 1; i <= n; i++) if (f[i] == top) {
            for (j = 1; j <= n; j++) if (map[i][j]) degree[j]--;
            degree[i] = -1;
        }
    }
    top--;
}

7. 网络预流和最大流
/***********************************************
网络中求最大流 Edmonds_Karp 最短增广路算法O(VE^2)
参数含义： n 代表网络中节点个数， 第1节点为源点， 第n节点为汇点
            net[][] 代表剩余网络， 0 表示无通路
            path[] 保存增广路径
            neck[]代表瓶颈， 保存增广路径最小容量
返回值： 最大流量
***********************************************/

const int NMAX = 210;
int net[NMAX][NMAX];
int path[NMAX], n;

int bfs() {
    queue<int> SQ;
    int neck[NMAX], i;
    memset(path, -1, sizeof(path));
    neck[1] = INT_MAX;
    SQ.push(1);
    while(!SQ.empty()) {
        int now = SQ.front();
        SQ.pop();
        if (now == n) break;
        for (i = 1; i <=n; ++i) if (net[now][i] > 0 && path[i] == -1) {
            path[i] = now;
            neck[i] = min(neck[now], net[now][i]);
            SQ.push(i);
        }
    }
    if (path[n] == -1) return -1;
    return neck[n];
}

int Edmonds_Karp() {
   int now, step;
   int max_flow = 0;
   while((step = bfs()) != -1) {
    max_flow += step;
    now = n;
    while(now != 1) {
        int pre = path[now];
        net[pre][now] -= step;
        net[now][pre] += step;
        now = pre;
    }
   }
   return max_flow;
}


/***********************************************
网络中求最大流HLPP高度标号预流推进算法O(V^2*E^0.5)
参数含义： n代表网络中节点的个数， 第0节点代表源点，第n节点为汇点
            net[][] 代表剩余网络， 0表示无通路
            earn[] 代表各点盈余
            high[] 代表各点的高度
返回值：    最大流量
***********************************************/
const int NMAX = 110;
int earn[NMAX], net[NMAX][NMAX], high[NMAX];
int n, m;
queue<int> SQ;
void push(int u, int v) {
    int ex = min(earn[u], net[u][v]);
    earn[u] -= ex;
    net[u][v] -= ex;
    earn[v] += ex;
    net[v][u] += ex;
}

void relable(int u) {
    int i, mmin = INT_MAX;
    for (i = 0; i <= n; i++) {
        if (net[u][i] > 0 && high[i] >= high[u]) {
            mmin = min(mmin, high[i]);
        }
    }
    high[u] = mmin + 1;
}

void discharge(int u) {
    int i, vn;
    while(earn[u] > 0) {
        vn = 0;
        for (i = 0; i <= n && earn[u] > 0; i++) if (net[u][i] > 0 && high[u] == high[i] + 1) {
            push(u, i);
            vn++;
            if (i != n) SQ.push(i);
        }
        if (vn == 0) relable(u);
    }
}

void init_preflow() {
    int i;
    memset(high, 0, sizeof(high));
    memset(earn, 0, sizeof(earn));
    while(!SQ.empty()) SQ.pop();
    high[0] = n + 1;
    for (i = 1; i <= n; i++) {
        if (net[0][i] > 0) {
            earn[i] = net[0][i];
            earn[0] -= net[0][i];
            net[i][0] = net[0][i];
            net[0][i] = 0;
            if (i != n) SQ.push(i);
        }
    }

}

int high_label_preflow_push() {
    int i, j;
    init_preflow();
    while(!SQ.empty()) {
        int overp = SQ.front();
        SQ.pop();
        discharge(overp);
    }
    return earn[n];

}


// 带gap优化的高标预流
const int N = 128;
const int INF = 1 << 28;
class Edge{
    public:
    int u, v, cuv, cvu, flow;
    Edge() {}
    Edge(int cu, int cv, int ccu, int ccv) : u(cu), v(cv), cuv(ccu), cvu(ccv), flow(0) {}
    int other(int p) const { return p == u ? v : u;}
    int cap(int p) const {return p == u ?cuv - flow : cvu - flow;}
    void addFlow(int p, int f) {flow += (p == u ? f : -f);}
};

class NodeList {
    private:
    int level, next[N], index[2 * N], v;
    public:
    void clear(int cv) {v = cv; level = -1; memset(index, -1, sizeof(index));}
    void insert(int n, int h) {next[n] = index[h]; index[h] = n; level >?=h;}
    int remove();
    bool empty() const {return level < 0;}
};

int NodeList::remove() {
    int r = index[level]; index[level] = next[index[level]];
    while(level >= 0 && index[level] == -1) level--;
    return r;
}

class NetWork {
    private:
        vector<Edge> eg;
        vector<Edge*> net[N];
        int v, s, t;
        NodeList list;
        int h[N], hn[2 * N], e[N], cur[N];
        void initNet();
        void initFlow();
        void initHeight();
        void push(int);
        void relabel(int);
        void discharge(int);
        void gapHeuristic(int);
    public:
        bool build();
        int maxFlow(int, int);
};

void NetWork::gapHeuristic(int k) {
    if (hn[k] != 0 || k >= v + 1) return ;
    for (int i = 0; i < v; i++) if(h[i] > k && h[i] <= v && i != s) {
        hn[h[i]]--; hn[v + 1]++; h[i] = v + 1;
    }
}

void NetWork::initNet() {
    for (int i = 0; i <v; i++) net[i].clear();
    for (int i = eg.size() - 1;i >= 0; i--) {
        net[eg[i].u].push_back(&eg[i]);
        net[eg[i].v].push_back(&eg[i]);
    }
}


void NetWork::initHeight() {
    memset(h, 0, sizeof(h)); memset(hn, 0, sizeof(hn));
    memset(e, 0, sizeof(e)); e[s] = INF;
    for (int i = 0; i < v; ++i) h[i] = v;
    queue<int> Q; Q.push(t); h[t] = 0;
    while(!Q.empty()) {
        int p = Q.front(); Q.pop();
        for (int i = net[p].size() - 1; i >= 0; i--) {
            int u = net[p][i]->other(p), ec = net[p][i]->cap(u);
            if (ec != 0 && h[u] == v && u != s) { h[u] = h[p] + 1; Q.push(u);}
        }
    }
    for (int i = 0; i < v; i++) hn[h[i]]++;
}

void NetWork::initFlow() {
    initNet(); initHeight();
    for (int i = 0; i < v; i++) cur[i] = net[i].size() - 1;
    list.clear(v);
    for(; cur[s] >= 0; cur[s]--) push(s);
}

void NetWork::push(int u) {
    Edge* te = net[u][cur[u]];
    int ex = min(te->cap(u), e[u]), p = te->other(u);
    if (e[p] == 0 && p != t) list.insert(p, h[p]);
    te->addFlow(u, ex); e[u] -= ex; e[p] += ex;
}

void NetWork::relabel(int u) {
    int mh = 2 * v, oh = h[u];
    for (int i = net[u].size() - 1; i >= 0; i--) {
        int p = net[u][i]->other(u);
        if (net[u][i]->cap(u) != 0) mh <?= h[p] + 1;
    }
    hn[h[u]]--; hn[mh]++; h[u] = mh; cur[u] = net[u].size() - 1;
    gapHeuristic(oh);
}

void NetWork::discharge(int u) {
    while(e[u] > 0)
        if (cur[u] < 0) relabel(u);
        else if (net[u][cur[u]]->cap(u) > 0 && h[u] == h[net[u][cur[u]]->other(u)] + 1) push(u);
        else cur[u]--;
}

bool NetWork::build() {
    int m, np, nc;
    int a, b, l, i;
    if (scanf("%d%d%d%d", &v, &np, &nc, &m) != 4) return false;
    v += 2; eg.clear();
    for (i = 0; i < m; i++) {
        scanf("\n(%d, %d)%d", &a, &b, &l);
        eg.push_back(Edge(a + 2, b + 2, l, 0));
    }
    for (i = 0; i < np; i++) {
        scanf("\n(%d)%d", &a, &l);
        eg.push_back(Edge(0, a + 2, l, 0));
    }

    for (i = 0; i < nc; i++) {
        scanf("\n(%d)%d", &a, &l);
        eg.push_back(Edge(a + 2, 1, l, 0));
    }
    return true;
}

int NetWork::maxFlow(int ss, int tt) {
    s = ss; t = tt; initFlow();
    while(!list.empty()) {
        int u = list.remove();
        discharge(u);
    }
    return e[t];

}

int main() {
    NetWork net;
    while(net.build()) printf("%d\n", net.maxFlow(0, 1));
    return 0;
}


/****************************************************
网络中求最大流Dinic算法O(V^2E)
适用于稠密图， 实际复杂度低于HLPP模板
参数含义： n 代表网络中节点数， 第0节点为源点， 第n节点为汇点
            net 代表网络， 似乎用前向星表示法存储边
            dis[] 代表从源点出发的距离标号
            path[] 代表模拟栈中的路径信息
            cur[] 代表模拟栈的现场保存
返回值：    最大流量
***************************************************/

const int NMAX = 21000;
const int MMAX = 250000<<1;
struct EDGE {
    int u, v, cap, flow;
    int next;
    EDGE(int _u = 0, int _v = 0, int _c = 0, int _f = 0)
    : u(_u), v(_v), cap(_c), flow(_f) {}
};

const int ENDFLAG = -1;
struct EDGELIST {
    int start[NMAX];
    int last[NMAX];
    int tot;
    EDGE arc[MMAX];
    void clear() {
        tot = ENDFLAG + 1;
        memset(last, ENDFLAG, sizeof(last));
    }
    void push_back(EDGE edge) {
        edge.next = ENDFLAG;
        arc[tot] = edge;
        if (last[edge.u] != ENDFLAG) arc[last[edge.u]].next = tot;
        else start[edge.u] = tot;

        last[edge.u] = tot;
        tot++;
    }

    // 创建双向弧
    void add_arc(EDGE edge) {
        push_back(edge);
        push_back(EDGE(edge.v, edge.u, edge.cap));
    }
}net;

int que[2][NMAX];
int qf[2], qe[2], qnow;

#define push_que(a) (que[qnow][qe[qnow]++] = (a))
#define pop_que2    (que[qnow^1][qf[qnow^1]++])
#define switch_que qnow ^= 1; \
                    qf[qnow] = qe[qnow] = 0;
#define empty_que2 (qf[qnow^1] >= qe[qnow^1])
#define size_que2 (qe[qnow^1] - qf[qnow^1])

int n, m;
int dis[NMAX];
int path[NMAX], deep;
int cur[NMAX];
bool bfs() {
    int i, j;
    memset(dis, -1, sizeof(dis));
    dis[0];
    qnow = 0;
    switch_que;
    push_que(0);
    switch_que;
    while(!empty_que2) {
        int l = size_que2;
        while(l--) {
            int u = pop_que2;
            for (i = net.start[u]; i != ENDFLAG; i = net.arc[i].next) {
                int v = net.arc[i].v;
                if ( dis[v] == -1 && net.arc[i].cap > net.arc[i].flow) {
                    push_que(v);
                    dis[v] = dis[u] + 1;
                    if (v == n) return true;
                }
            }
        }
        switch_que;
    }
    return false;
}

int Dinic() {
    int i, j;
    int u;
    int maxflow = 0;
    while(bfs()) {
        memcpy(cur, net.start, sizeof(cur));
        for (deep = u = 0; true;) {
            if (u == n) {
                int neck = INT_MAX, pos;
                for (i = 0; i < deep; i++) {
                    int res = net.arc[path[i]].cap - net.arc[path[i]].flow;
                    if (res < neck) {
                        neck = res;
                        pos = i;
                    }
                }
                maxflow += neck;
                for (i = 0; i < deep; i++) {
                    net.arc[path[i]].flow += neck;
                    net.arc[path[i]^1].flow -= neck;
                }
                deep = pos;
                u = net.arc[path[deep]].u;
            }
            for (i = cur[u]; i != ENDFLAG; i = net.arc[i].next) {
                if (net.arc[i].cap > net.arc[i].flow &&
                    dis[u] + 1 == dis[net.arc[i].v]) break;
            }
            cur[u] = i;
            if (i != ENDFLAG) {
                path[deep++] = i;
                u = net.arc[i].v;
            } else {
                if (deep == 0) break;
                dis[u] = -1;
                u = net.arc[path[--deep]].u;
            }
        }
    }
    return maxflow;
}

8. 网络最小费用最大流
/*********************************************
网络中最小费用最大流 O(V*E^2)
参数含义： n 代表网络中的节点数， 第0节点为源点， 第n节点为汇点
            net[][] 代表剩余网络
            cost[][] 代表单位费用
            path[] 保存增广路径
            ecost[] 源点到各点的最短路

算法： 初始最小费用和最大流均为0， 寻找单位费用最短路
在最短路中求出最大流， 即为增广路， 再修改剩余网络，直到无可增广路为止。
返回值： 最小费用， 最大流量
****************************************/

const int NMAX = 210;
int net[NMAX][NMAX], cost[NMAX][NMAX];
int path[NMAX], ecost[NMAX];
int n;
bool bellman_ford() {
    int i, j;
    memset(path, -1, sizeof(path));
    fill(ecost, ecost + NMAX, INT_MAX);
    ecost[0] = 0;

    bool flag = true;
    while(flag) {
        flag = false;
        for (i = 0; i <= n; i++) {
            if (ecost[i] == INT_MAX) continue;
            for (j = 0;j <= n; j++) {
                if (net[i][j] > 0 && ecost[i] + cost[i][j] < ecost[j]) {
                    flag = true;
                    ecost[j] = ecost[i] + cost[i][j];
                    path[j] = i;
                }
            }
        }
    }
    return ecost[n] != INT_MAX;
}

int min_cost_max_flow() {
    int i, j;
    int mincost = 0, maxflow = 0;
    while(bellman_ford()) {
        int now = n;
        int neck = INT_MAX;
        while(now != 0) {
            int pre = path[now];
            neck = min(neck, net[pre][now]);
            now = pre;
        }
        maxflow += neck;
        now = n;
        while(now != 0) {
            int pre = path[now];
            net[pre][now] -= neck;
            net[now][pre] += neck;
            cost[now][pre] = -cost[pre][now];
            mincost += cost[pre][now] * neck;
            now = pre;
        }
    }
    return mincost;
}
/*****************************************************
网络中最小费用最大流O(V*E^2) 邻接表SPFA实现
参数含义： n 代表网络中节点的个数， 第s节点为源点， 第t节点为汇点
            net代表剩余网络
            path 保存增广路径
            ecost 保存源点到各点的最短路

返回值：    最小费用， 最大流量
*****************************************************/
// POJ 3422
const int NMAX = 5100;
const int MMAX = 30000;
const int INF = 0x7f7f7f7f;

int path[NMAX], ecost[NMAX];
int n;
int s, t;

struct EDGE {
    int u, v, cap, cost, flow;
    int next;
    EDGE(int _u=0, int _v = 0, int _c = 0, int _ct = 0, int _f = 0):
        u(_u), v(_v), cap(_c), flow(_f), cost(_ct) {}
};

const int ENDFLAG = -1;
struct EDGELIST {
    int start[NMAX];
    int last[NMAX];
    int tot;
    EDGE arc[MMAX];
    void clear() {
        tot = ENDFLAG + 1;
        memset(last, ENDFLAG, sizeof(last));
    }
    void push_back(EDGE edge) {
        edge.next = ENDFLAG;
        arc[tot] = edge;
        if (last[edge.u] != ENDFLAG) arc[last[edge.u]].next = tot;
        else start[edge.u] = tot;
        last[edge.u] = tot;
        tot++;
    }
    // 创建双向弧
    void add_arc(EDGE edge) {
        push_back(edge);
        push_back(EDGE(edge.v, edge.u, 0, INF));
    }
}net;

int que[2][NMAX];
int qf[2], qe[2], qnow;

#define push_que(a) (que[qnow][qe[qnow]++] = (a))
#define pop_que2    (que[qnow^1][qf[qnow^1]++])
#define switch_que qnow ^= 1; \
                    qf[qnow] = qe[qnow] = 0;
#define empty_que2 (qf[qnow^1] >= qe[qnow^1] )
#define size_que2 (qe[qnow^1] - qf[qnow^1])

bool SPFA() {
    int i, j;
    bitset<NMAX> vis;
    memset(ecost, 0x7f, sizeof(ecost));
    memset(path, -1, sizeof(path));

    bool flag = true;
    qnow = 1;
    switch_que;
    push_que(s);
    vis[s] = 1;
    ecost[s] = 0;
    for (j = 0; j < n && flag; j++) {
        flag = false;
        switch_que;
        int l = size_que2;
        while(l--) {
            int now = pop_que2;
            vis[now] = 0;
            for (i = net.start[now]; i != ENDFLAG; i = net.arc[i].next) {
                EDGE ed = net.arc[i];
                if (ed.cap > ed.flow && ecost[ed.v] > ecost[now] + ed.cost) {
                    flag = true;
                    ecost[ed.v] = ecost[now] +ed.cost;
                    path[ed.v] = i;
                    if (!vis[ed.v]) {
                        vis[ed.v] = 1;
                        push_que(ed.v);
                    }
                }
            }
        }
    }
    return ecost[t] != INF;
}

int min_cost_max_flow() {
    int i, j;
    int mincost = 0, maxflow = 0;
    while(SPFA()) {
        int pre = path[t];
        int neck = INT_MAX;
        while(pre != -1) {
            int res = net.arc[pre].cap - net.arc[pre].flow;
            neck = min(neck, res);
            pre = path[net.arc[pre].u];
        }
        maxflow += neck;
        mincost += ecost[t] * neck;
        pre = path[t];
        while(pre != -1) {
            net.arc[pre].flow += neck;
            net.arc[pre^1].flow -= neck;
            net.arc[pre^1].cost = -net.arc[pre].cost;
            pre = path[net.arc[pre].u];
        }
    }
    return mincost;
}

9. 网络最大流（高度标号预流推进）
/**************
函数接口： int Relabel_To_Front(int s, int d) O(V^2*sqrt(E))
参数含义：   s为源点， d为汇点
返回值： 网络最大流

调用函数前的初始化工作： ver 置为网络中节点的个数， c[i][j]代表节点i到节点j的流量，
vl[i] 存放i与相邻的所有节点。
其他全局变量都初始化为0
*/

const int VEX = 405; // 网络中顶点的个数
const int HMAX = 810; // 最大高度的定义，只要大于顶点的2倍就可以了
int f[VEX][VEX]; // 流量
int c[VEX][VEX]; // 边最大容量
int h[VEX]; // 节点高度
int e[VEX]; // 节点容量
int ver; // 节点数目

vector<int> vl[VEX]; // 邻接表， vl[i]存放与i相邻的节点。

void Push(int u, int v) {// 流推进， 由节点u推向节点v
    int cf = c[u][v] - f[u][v]; // u, v边的容量
    int d = e[u] < cf ? e[u] : cf;
    f[u][v] += d;
    f[v][u] = -f[u][v];
    e[u] -= d;
    e[v] += d;
}

void Relabel(int u) {// 对u重新标号
    int i, t, cf;
    int hmin = HMAX;
    for (i = 0; i < vl[u].size(); i++) { // 寻找相邻最低点
        t = vl[u][i];
        cf = c[u][t] - f[u][t];
        if (cf > 0 && h[u] <= h[t] && h[t] < hmin) hmin = h[t];
    }
    h[u] = hmin + 1;
}

void Init_Preflow(int s) { // 初始化网络流， s为源点
    int i;
    int u;
    h[s] = ver; //初始化高度
    for(i = 0; i < vl[s].size(); i++) {
        u = vl[s][i];
        f[s][u] = c[s][u];
        f[u][s] = -c[s][u];
        e[u] = c[s][u];
        e[s] -= c[s][u];
    }
}

void Discharge(int u) {
    int i = 0;
    int cf, v;
    if (vl[u].size() == 0) return ;
    while(e[u] > 0) {
        if (i < vl[u].size()) {
            v = vl[u][i];
            cf = c[u][v] - f[u][v];
        }
        if (i >= vl[u].size()) {
            Relabel(u);
            i = 0;
        } else if (cf > 0 && h[u] == h[v] + 1)
            Push(u, v);
        else i++;
    }
}

int Relabel_To_Front(int s, int d) {// s 为源点， d为汇点
    int u, i, old_h;
    list<int> l;
    list<int>::iterator iter;

    Init_Preflow(s);

    iter = l.begin();
    for (i = 0; i < ver; i++) {
        if (i != s && i != d) l.insert(iter, i);
    }
    iter = l.begin();
    while(iter != l.end()) {
       u = *iter;
       old_h = h[u];
       Discharge(u);
       if (h[u] > old_h) {
        l.erase(iter);
        l.insert(l.begin(), u);
        iter = l.begin();
       }
       iter++;
    }
    return e[ver - 1];
}

10. 最大团
/**********************************************
最大独立集，最大团
PKU 1419 Graph Coloring
团： 指G的一个完全子图， 该子图不包含在任何其他的完全子图当中
最大独立集： 补图的最大团
最大团： 指其中包含顶点最多的团
**********************************************/
#include <cstdio>
#include <string>
#define NMAX 110
bool path[NMAX][NMAX];
int n, mmax;
int dp[NMAX];
bool v[NMAX];
int seq[NMAX], seq_pos;
// seq 记录最大团集合
bool dfs(int pos, int size) {
    int i, j, unvis;
    bool tv[NMAX];
    unvis = 0;
    for (i = pos; i < n; i++) if (!v[i]) {
        unvis++;
    }
    if (unvis == 0) {
        if (size > mmax) {
            mmax = size;
            seq_pos = 0;
            seq[seq_pos++] = pos + 1;
            return true;
        }
        return false;
    }
    for (i = pos; i < n && unvis > 0; i++) {
        if (!v[i]) {
            if (unvis + size <= mmax || dp[i] + size <= mmax) return false;
            v[i] = true;
            unvis--;
            memcpy(tv, v, sizeof(v));
            for (j = 0; j < n; j++) {
                if (!path[i][j]) return true;
            }
            if (dfs(i, size + 1)) {
                seq[seq_pos++] = pos + 1;
                return true;
            }
            memcpy(v, tv, sizeof(v));
        }
    }
    return false;
}

int max_clique() {
    int i, j;
    mmax = 0;
    for (i = 0; i < n; i++) path[i][i] = false;
    for (i = n - 1; i >= 0; i--) {
        for (j = 0; j < n; j++) v[j] = !path[i][j];
        dfs(i, 1);
        dp[i] = mmax;
    }
    return mmax;
}

int main() {
    int i, j, x,y, e;
    int m, tn;
    scanf("%d", &m);
    while(m--) {
        scanf("%d %d", &n, &e);
        memset(path, 0, sizeof(path));
        for (i = 0; i < e; i++) {
            scanf("%d %d", &x, &y);
            x--; y--;
            path[x][y] = path[y][x] = true;
        }
        // max independent set in original graph
        // max clique in inverse graph
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                path[i][j] = !path[i][j];
            }
        }
        memset(dp, 0, sizeof(dp));
        printf("%d\n", max_clique());

        printf("%d", seq[0]);
        for (i = 1; i < seq_pos; i++) printf(" %d", seq[i]);
        printf("\n");
    }
}

11. 最大二分图匹配（匈牙利算法）
/************************************
二分图： 指所有的顶点分成集合M和N， M或者N中任意两个在同一集合里的点互不相连
匹配： 一组边顶点分别在两个集合中， 并且任意两条边都没有相同顶点
最大匹配： 所能得到的最大的边的个数
************************************/
#include <cstdio>
#include <memory>
#include <vector>
using namespace std;
const int Max = 1100;
vector< vector<int> > Bmap;
int n, m, k, nm;
int mark[Max];
bool flag[Max];

bool dfs(int pos) {
    int i, pre, tp;
    for (i = 0; i < Bmap[pos].size(); i++) {
        tp = Bmap[pos][i];
        if (!flag[tp]) {
            flag[tp] = true;
            pre = mark[tp];
            mark[tp] = pos;
            if (pre == -1 || dfs(pre)) return true;
            mark[tp] = pre;
        }
    }
    return false;
}

inline int Max_Match() {
    int mmax = 0, i;
    for (i = 1; i <= m; i++) {
        memset(flag, 0, sizeof(flag));
        if (dfs(i)) mmax++;
    }
    return mmax;
}

int main() {
    int i, j, id, id2;
    while(scanf("%d", &k) == 1 && k) {
        scanf("%d%d", &m, &n);
        nm = n + m;
        Bmap.clear(); Bmap.resize(nm + 10);
        memset(mark, -1, sizeof(mark));
        for (j = 0; j < k; j++) {
            scanf("%d %d", &id, &id2);
            id2 += m;
            Bmap[id].push_back(id2);
        }
        printf("%d\n", Max_Match());
    }
}

// 二分匹配HopcroftKarp 算法 O(sqrt(V)*E)
// 贪心一个初始匹配可以加速
#include <iostream>
#include <queue>
using namespace std;

const int MAXN = 3002;
const int INF = 1<<30;

struct node {
    int x, y;
}G[MAXN], U[MAXN];

int n, m, t, nx, ny, dis;
int x[MAXN], y[MAXN], vs[MAXN];
int Isf[MAXN];
bool adj[MAXN][MAXN];
int ds[MAXN], dt[MAXN];

void input() {
    scanf("%d %d", &t, &m);
    for (int i = 0; i < m; i++) {
        scanf("%d %d %d", &G[i].x, &G[i].y, &vs[i]);
        vs[i] *= vs[i];
    }
    scanf("%d", &n);
    for (int i = 0; i < n; i++) scanf("%d %d", &U[i].x, &U[i].y);
    memset(adj, 0, sizeof(adj));
    if (m < n) {
        nx = m, ny = n;
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++) {
                int a = G[i].x - U[j].x,  b = G[i].y - U[j].y;
                if ((a * a + b * b) < vs[i] * t * t) adj[i][j] = 1;
            }
    } else {
        nx = n, ny = m;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) {
                int a = G[j].x - U[i].x, b = G[j].y - U[i].y;
                if ((a * a + b * b) <= vs[j] * t * t) adj[i][j] = 1;
            }
    }
}

bool Search() {
    memset(ds, -1, sizeof(ds));
    memset(dt, -1, sizeof(dt));
    queue<int> Q; dis = INF;
    for (int i = 0; i < nx; i++) if (x[i] == -1) {
        Q.push(i);
        ds[i] = 0;
    }
    while(!Q.empty()) {
        int u = Q.front(); Q.pop();
        if (ds[u] > dis) break;
        for (int v = 0; v < ny; v++) {
            if (adj[u][v]) {
                if (dt[v] != -1) continue;
                dt[v] = ds[u] + 1;
                if (y[v] == -1) dis = dt[v];
                else {
                    ds[y[v]] = dt[v] + 1;
                    Q.push(y[v]);
                }
            }
        }
    }
    return (dis != INF);
}

bool DFS(int u) {
    for (int v = 0; v < ny; v++) {
        if (Isf[v] && adj[u][v] && dt[v] == ds[u] + 1) {
            Isf[v] = true;
            if (y[v] != -1 && dt[v] == dis) continue;
            if (y[v] == -1 || DFS(y[v])) {
                y[v] = u; x[u] = v;
                return true;
            }
        }
    }
    return false;
}

int HopcroftKarp() {
    int cnt = 0;
    for (int i = 0; i < nx; i++) x[i] = i;
    for (int i = 0; i < ny; i++) y[i] = -1;
    while(Search()) {
        memset(Isf, 0, sizeof(Isf));
        for (int i = 0; i < nx; i++) if (dis[i] == 0 && DFS(i)) cnt++;
    }
    return cnt;
}

int main() {
    int test;
    scanf("%d", &test);
    for (int k = 1; k <= test; k++) {
        input();
        printf("Scenario #%d:\n%d\n\n", k, HopcroftKarp());
    }
    return 0;
}

12. 带权二分图最优匹配（KM算法）
// 二分图带权匹配O(N^3)
const int MAXN = 509;
const int INF = 0x1fffffff;

int bpCostMatch(int c[][MAXN], int nx, int ny) {
    static int lx[MAXN], ly[MAXN], slack[MAXN];
    static int open[MAXN], prev[MAXN], pnt[MAXN], x[MAXN], y[MAXN];
    int i, j, k, s, head, tail;
    int d, ans = 0;

    if (nx > ny) ny = nx;
    for (i = 0; i < nx; i++) lx[i] = -INF;
    for (i = 0; i < ny; i++) ly[i] = 0;
    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++) if ((lx[i] - c[i][j]) < 0) {
            lx[i] = c[i][j];
        }
    memset(x, -1, sizeof(x)); memset(y, -1, sizeof(y));
    for (i = 0; i < nx; i++) {
        memset(prev, -1, sizeof(prev));
        for (j = 0; j < ny; j++) slack[j] = INF;
        open[0] = i; head = 0; tail = 1;
        while(x[i] < 0) {
            for (; head < tail && x[i] < 0; head++)
                for (s = open[head], j = 0; j < ny && x[i] < 0; j++) if (prev[j] < 0) {
                    if ((d = lx[s] + ly[j] - c[s][j]) > 0) {
                        if ((slack[j] - d) > 0) {
                            slack[j] = d; pnt[j] = s;
                        }
                        continue;
                    }
                    open[tail++] = y[j]; prev[j] = s;
                    if (y[j] >= 0) continue;
                    while(j >= 0) {
                        s = prev[j]; y[j] = s; k = x[s]; x[s] = j; j = k;
                    }
                }
                if (x[i] >= 0) break;
                for (d = INF, j = 0; j < ny; j++) if (prev[j] < 0 && (d - slack[j]) > 0) d = slack[j];
                for (j = 0; j < tail; j++) lx[open[j]] -= d;
                for (j = 0; j < ny; j++) if (prev[j] >= 0) {
                    ly[j] += d;
                } else if(slack[j] < INF) {
                    slack[j] -= d;
                }
                for (j = 0; j < ny; j++)
                    if (prev[j] < 0 && slack[j] == 0) {
                        open[tail++] = y[j]; prev[j] = pnt[j];
                        if (y[j] >= 0) continue;
                        while(j >= 0) {
                            s = prev[j]; y[j] = s; k = x[s]; x[s] = j;j = k;
                        }
                        break;
                    }

        }
    }
    for (i = 0; i < nx; i++)
        if (c[i][x[i]] > -INF) {
            if (c[i][x[i]] < 0)
                return -1;
            ans += c[i][x[i]];
        } else return -1;
        return ans;
}

int N, M, E;
int c[MAXN][MAXN];

int cas;
int main() {
    int i, j, a, b, w, ans;
    while(scanf("%d%d%d", &N, &M, &E) != EOF) {
        for (i = 0; i < N; i++)
            for (j = 0; j < M; j++)
                c[i][j] = -INF;
        for (i = 0; i < E; i++) {
            scanf("%d%d%d", &a, &b, &w);
            if (w < 0) continue;
            if (c[a][b] < w) c[a][b] = w;
        }
        if (N > M) ans = -1;
        else ans = bpCostMatch(c, N, M);
        printf("Case %d:%d\n", ++cas, ans);
    }
    return 0;
}


/*
wywcgs 的KM O(n^3)
只需要把Graph里的n（顶点数） 和 edgep[x][y](边权)赋值， 第一维为x点，
第二维为y点。
然后调用KMMatch()函数即可，返回值为最大权完美匹配。
最小权匹配可将每条边取反， 然后类似求最大权匹配即可。
匹配信息保存在xmate[] 和 ymate[] 中。 其中xmate[i] 为x[i]的匹配点，
ymate[i]为y[i]的匹配点。
*/
#include <cstdio>
#include <queue>
#include <algorithm>
using namespace std;

const int N = 310;
const int INF = 1<<28;

class Graph {
    private:
        bool xckd[N], yckd[N];
        int n, edge[N][N], xmate[N], ymate[N];
        int lx[N], ly[N], slack[N], prev[N];
        queue<int> Q;
        bool bfs();
        void agument(int);
    public:
        bool make();
        int KMMatch();
};

bool Graph::bfs() {
    while(!Q.empty()) {
        int p = Q.front(), u = p>>1; Q.pop();
        if (p&1) {
            if (ymate[u] == -1) {agument(u); return true;}
            else {xckd[ymate[u]] = true; Q.push(ymate[u]<<1);}
        } else {
            for (int i = 0; i < n; i++)
                if (yckd[i]) continue;
                else if (lx[u] + ly[i] != edge[u][i]) {
                    int ex = lx[u] + ly[i] - edge[u][i];
                    if (slack[i] > ex) {slack[i] = ex; prev[i] = u;}
                } else {
                    yckd[i] = true; prev[i] = u;
                    Q.push((i<<1) | 1);
                }
        }
    }
    return false;
}

void Graph::agument(int u) {
    while(u != -1) {
        int pv = xmate[prev[u]];
        ymate[u] = prev[u]; xmate[prev[u]] = u;
        u = pv;
    }
}

int Graph::KMMatch() {
    memset(ly, 0, sizeof(ly));
    for (int i = 0; i < n; i++) {
        lx[i] = -INF;
        for (int j = 0; j < n; j++) lx[i] >?= edge[i][j];
    }
    memset(xmate, -1, sizeof(xmate)); memset(ymate, -1, sizeof(ymate));
    bool agu = true;
    for (int mn = 0; mn < n; mn++) {
        if (agu) {
            memset(xckd, false, sizeof(xckd));
            memset(yckd, false, sizeof(yckd));
            for (int i = 0; i < n; i++) slack[i] = INF;
            while(!Q.empty()) Q.pop();
            xckd[mn] = true; Q.push(mn<<1);
        }
        if (bfs()) {agu = true; continue;}
        int ex = INF; mn--; agu = false;
        for (int i = 0; i < n; i++) if (!yckd[i]) ex<?=slack[i];
        for (int i = 0; i <n; i++) {
            if(xckd[i]) lx[i] -= ex;
            if (yckd[i]) ly[i] += ex;
            slack[i] -= ex;
        }
        for (int i = 0; i < n; i++)
            if (!yckd[i] && slack[i] == 0)  { yckd[i] = true; Q.push((i<<1) | 1);}
    }
    int cost = 0;
    for (int i = 0; i < n; i++) cost += edge[i][xmate[i]];
    return cost;
}

bool Graph::make() {
    int i, j;
    while(scanf("%d", &n) == 1) {
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                scanf("%d", &edge[i][j]);
        return true;
    }
    return false;
}

int main() {
    Graph g;
    while(g.make()) printf("%d\n", g.KMMatch());
    return 0;
}

13. 强连通分量(Kosaraju算法)
/*
有向图的强连通分量Kosaraju算法O(E)
参数含义： 使用邻接表来保存图
            path原图， npath逆图
            scc强连通的个数
            id[x] = y 表示第x个顶点属于y强连通
*/
#define NMAX 11000
vector< vector<int> > path;
vector< vector<int> > npath;
int n, m, scc;
int order[NMAX], order_pos, id[NMAX];
bool vis[NMAX];

void dfs(int pos) {
    int i, j, l;
    vis[pos] = true;
    l = path[pos].size();
    for (i = 0; k < l; i++) {
         j = path[pos][i];
         if (!vis[j]) dfs(j);
    }
    order[order_pos++] = pos;
}

void ndfs() {
    int i, j, l;
    vis[pos] = true;
    id[pos] = scc;
    l = npath[pos].size();
    for (i = 0; i < l; i++) {
        j = npath[pos][i];
        if (!vis[j]) ndfs(j);
    }
}

void Kosaraju() {
    int i, j;
    // dfs in original graph
    memset(vis, 0, sizeof(vis));
    order_pos = 0;
    for (i = 1; i <= n; i++) {
        if (!vis[i]) dfs(i);
    }
    // dfs in inverse graph
    memset(vis, 0, sizeof(vis));
    memset(id, 0, sizeof(id));
    scc = 1;
    for (i = order_pos - 1; i >= 0; i--) if (!vis[order[i]]) {
        ndfs(order[i]);
        scc++;
    }
    scc--;
}

14. 强连通分量（Gabow算法）
/*
有向图的强连通分量Gabow算法O(E)
参数含义： 使用邻接表来保存图
            path原图
            scc强连通个数
            id[x] = y 表示第x个顶点属于y强连通
*/
#define NMAX 11000
vector< vector<int> > path;
int n, m, scc, step;
int order[NMAX], order_pos, id[NMAX];
int order2[NMAX], order2_pos;
int vis[NMAX];

void dfs(int pos) {
    int i, j, next, l, pre;
    vis[pos] = step++;
    order[order_pos++] = pos;
    order2[order2_pos++] = pos;
    l = path[pos].size();
    for (i = 0; i < l; i++) {
        next = path[pos][i];
        if (vis[next] == 0) dfs(next);
        else if (id[next] == 0) { // have a circle and belongs to nothing
            while(vis[order2[order2_pos - 1]] > vis[next]) order2_pos--;
        }
    }
    if (order2[order2_pos - 1] == pos) order2_pos--;
    else return ;
    do {
        pre = order[order_pos - 1];
        id[pre] = scc;
        order_pos--;
    } while(pre != pos);
    scc++;
}

void Gabow() {
    int i, j, l;
    // dfs in original graph
    memset(id, 0, sizeof(id));
    memset(vis, 0, sizeof(vis));
    scc = step = 1;
    order_pos = order2_pos = 0;
    for (i = 1; i <= n; i++) if (vis[i] == 0) dfs(i);
    scc--;
}

15. 无向图割边割点和双连通分量
#define mclear(x) memset((x), 0, size((x)))
const int MAX = 5100;
int n, m, deep;
vector<int> path[MAX];
int vis[MAX], low[MAX];
vector<int> cutpoint; // 割点
vector< pair<int, int> > bridge; // 割边
int nbcc; //双连通分量数
stack< pair<int, int> > order;
vector<int> bcc[MAX]; // 双连通分量

void dfs(int pos, int father) {
    int i, j, total = 0;
    bool cut = false;
    int reback = 0; // 处理平行边
    vis[pos] = low[pos] = deep++;
    int ls = path[pos].size();
    for (j = 0; j < ls; j++) {
        i = path[pos][j];
        if (i == father) reback++;
        if (vis[i] == 0) {
            pair<int, int> e(pos, i);
            order.push(e);
            dfs(i, pos);
            if (low[i] >= vis[pos]) {
                nbcc++;
                bcc[nbcc].clear();
                pair<int, int> r;
                do {
                    r = order.top();
                    order.pop();
                    bcc[nbcc].push_back(r.second);
                } while(e != r);
                bcc[nbcc].push_back(r.first);
            }
            total++;
            low[pos] = min(low[i], low[pos]);
            if ((vis[pos] == 1 && total > 1) ||
                (vis[pos] != 1 && low[i] >= vis[pos])) cut = true;
            if (low[i] > vis[pos]) bridge.push_back(e);
        } else if (i != father) {
            low[pos] = min(vis[i], low[pos]);
        }
    }
    if (reback > 1) low[pos] = min(low[pos], vis[father]);
    if (cut) cutpoint.push_back(pos);
}

void find_cut() {
    int i;
    mclear(vis); mclear(low);
    cutpoint.clear(); bridge.clear();
    nbcc = 0;
    while(!order.empty()) order.pop();
    for (i = 1; i <= n; i++) if (vis[i] == 0) {
        deep = 1;
        dfs(i, i);
    }
}

/********************************************************
图的DFS信息构建 By oyjpArt
g矩阵： g[i][j] -> 0 : 无边
                -> 1 : 可重复访问边
                -> -1: 非可重复访问边
说明： 以为在无向图中u->v 访问之后就不能再从v->u访问了
故{u, v} 访问了之后 {v, u}要置为-1
如果是有向图则没有这个规则
gc矩阵：gc[i][j] -> 0 : 无边
                -> 1 : 树枝边
                -> 2 : 反向边
                -> 3 : 正向边
                -> 4 : 交叉边
d数组：    顶点的开始访问时间表
f数组：    顶点的结束访问时间表
c数组：    顶点颜色表 0 白色， -1 灰色， 1 黑色
p数组：    顶点的前驱表
l数组：    顶点的L值（最顶层的祖先层数）
b数组：    顶点的度数表

关于标号函数LOW()
LOW(U) 代表的是与U以及U的子孙直接相连的节点的最高辈分（深度）
d[U]  U首次被访问时
LOW[U] = min(LOW[U], d[W]); 访问边{U, W}
min(LOW[U], LOW[S]) U的儿子S的关联边全部被访问是
********************************************************/

const int maxn = 100;
int n, g[maxn][maxn], gc[maxn][maxn];
int d[maxn], f[maxn], l[maxn], p[maxn], c[maxn], b[maxn];
int time;

void dfs_visit(int u) { // 递归搜索以U为根的深度优先树
    int v;
    c[u] = -1;
    time++; d[u] = time; l[u] = time;
    for (v = 1; v <= n; v++) if (g[u][v] > 0) {
        if (c[v] == 0) {
            g[v][u] = -1;
            gc[u][v] = 1;
            b[u]++;
            p[v] = u;
            dist_visit(v);
            if (l[v] <l[u])
                l[u] = l[v];
            g[v][u] = 1;
        } else {
            if (c[v] < 0) {
                if(l[v] < l[u])
                    l[u] = l[v];
                gc[u][v] = 2;
            } else {
                if (d[v] > d[u])
                    gc[u][v] = 3;
                else gc[u][v] = 4;
            }
        }

    }
    c[u] = 1;
    time++; f[u] = time;

}

void dfs() {
    int u;
    memset(gc, 0, sizeof(gc));
    memset(c, 0, sizeof(c));
    memset(b, 0, sizeof(b));
    time = 0;
    for (u = 1; u <= n; u++) if (c[u] == 0) {
        p[u] = 0;
        dfs_visit(u);
    }
}

16. 最小树形图O(N^3)
/*
最小树形图 O(N^3)
参数含义： 使用邻接表来保存图， 邻接表O(VE)
            path 原图
            pre 保存最小入弧的权
            del 表示被缩去的点
            fpre 保存最小树形图的逆路径
例题： TJU 2248
*/

const int NMAX = 110;
const int INF = 0x7f7f7f7f;
int n;
int path[NMAX][NMAX];
int pre[NMAX];
bool vis[NMAX], del[NMAX];
int min_cost;
int fold[NMAX], fpre[NMAX];
void dfs(int pos) {
    int i;
    vis[pos] = true;
    for (i = 1; i <= n; i++) if (path[pos][i] != INF && !vis[i]) dfs(i);
}

bool is_connect(int root) {
    int i;
    memset(vis, 0, sizeof(vis));
    dfs(root);
    for (i = 1; i <= n; i++) if (!vis[i]) return false;
    return true;
}

//O(N^3)
bool min_tree_graph(int root) {
    int i, j, k;
    // make sure every node(except root) have in-arc
    if (!is_connect(root)) return false;
    memset(del, 0, sizeof(del));
    min_cost = 0;
    for (i = 0; i <= n; i++) fold[i] = fpre[i] = i;
    while(true) {
        for (i = 1; i <= n; i++) {
            if (del[i] || i == root) continue;
            pre[i] = i;
            path[i][i] = INF; // delete self-cycle
            for (j = 1; j <= n; j++) {
                if (del[j]) continue;
                if (path[j][i] < path[pre[i]][i]) pre[i] = fpre[fold[i]] = j;
            }
        }
        for (i = 1; i <= n; i++) {
            if (del[i] || i == root) continue;
            j = i;
            memset(vis, 0, sizeof(vis));
            while(!vis[j] && j != root) {
                vis[j] = true;
                j = pre[j];
            }
            if (j == root) continue;
            i = j;
            min_cost += path[pre[i]][i];
            for (j = pre[i]; j != i; j = pre[j]) {
                del[j] = true;
                min_cost += path[pre[j]][j];
            }
            for (j = 1; j <= n; j++) {
                if (del[j]) continue;
                if (path[j][i] != INF) path[j][i] -= path[pre[i]][i];
            }
            for (j = pre[i]; j != i; j = pre[j]) {
                for (k = 1; k <= n; k++) {
                    if (del[k]) continue;
                    path[i][k] = min(path[i][k], path[j][k]);
                    if (path[k][j] != INF && path[k][i] > path[k][j] - path[pre[j]][j]) {
                        path[k][i] = path[k][j] - path[pre[j]][j];
                        fold[i] = j;
                        fpre[i] = j;
                    }
                }
            }
            break;
        }
        if (i > n) {
            for (i = 1; i <= n; i++) {
                if (del[i] || i == root) continue;
                min_cost += path[pre[i]][i];
            }
            break;
        }
    }
    return true;
}


// print path in min tree Graph
void print_mtg(int root) {
    int i, total = n;
    memset(vis, 0, sizeof(vis));
    for (i = 1; i <= n; i++) vis[fpre[i]] = true;
    for (i = 1; i <= n; i++) {
        if (!vis[i]) {
            int pos = i;
            while(pos != root) {
                printf("%d <- ", pos);
                pos = fpre[pos];
            }
            printf("%d\n", root);
        }
    }
}

int main() {
    int i, m;
    while(scanf("%d %d", &n, &m), !(n == 0 && m == 0)) {
        memset(path, ox7f, sizeof(path));
        while(m--) {
            int x, y, z;
            scanf("%d %d %d", &x, &y, &z);
            path[x][y] = min(path[x][y], z);
        }
        if (!min_tree_graph(1)) puts("impossilbe");
        else printf("%d\n", min_cost);
    }
}

17. 最小树形图 O(VE)
const int NMAX = 1500;
const int INF = 0x7f7f7f7f;
struct LINKT {
    int ls;
    int adj[NMAX];
    void clear() {ls = 0;}
    int operator[] (const int pos) {return adj[pos];}
    int size() {return ls;}
    void push_back(const int pos) { adj[ls++] = pos;}
};

int n;
int path[NMAX][NMAX];
LINKT epath[NMAX], nepath[NMAX];
int pre[NMAX];
bool vis[NMAX], del[NMAX];
int min_cost;
int fold[NMAX], fpre[NMAX];
void dfs(int pos) {
    int i ;
    vis[pos] = true;
    for (i = 0; i < epath[pos].ls; i++) if (!vis[epath[pos].adj[i]]) dfs(epath[pos].adj[i]);
}

bool is_connect(int root) {
    int i;
    memset(vis, 0, sizeof(vis));
    dfs(root);
    for (i = 1; i <= n; i++) if (!vis[i]) return false;
    return true;
}

// O(VE)
bool min_tree_graph(int root) {
    int i, j, k;
    // make sure every node(except root) have in-arc
    if (!is_connect(root)) return false;
    memset(del, 0, sizeof(del));
    min_cost = 0;
    for (i = 0; i <= n; i++) fold[i] = fpre[i] = i;
    while(true) {
        for (i = 1; i <= n; i++) {
            if (del[i] || i == root) continue;
            pre[i] = i;
            path[i][i] = INF;
            for (j = 0; j < nepath[i].ls; j++) {
                int t = nepath[i].adj[j];
                if (del[t]) continue;
                if (path[t][i] < path[pre[i]][i]) pre[i] = fpre[fold[i]] = t;
            }
        }
        for (i = 1; i <= n; i++) {
            if (del[i] || i == root) continue;
            j = i;
            memset(vis, 0, sizeof(vis));
            while(!vis[j] && j != root) {
                vis[j] = true;
                j = pre[j];
            }
            if (j == root) continue;
            i = j;
            min_cost += path[pre[i]][i];
            for (j = pre[i]; j != i; j = pre[j]) {
                del[j] = true;
                min_cost += path[pre[j]][j];
            }

            for (j = 0; j < nepath[i].ls; j++) {
                int t = nepath[i].adj[j];
                if (del[t]) continue;
                path[t][i] -= path[pre[i]][i];
            }

            for (j = pre[i]; j != i; j = pre[j]) {
                for (k = 0; k < epath[j].ls; k++) {
                    int t = epath[j].adj[k];
                    if (del[t]) continue;
                    if (path[i][t] == INF) {
                        epath[i].push_back(t);
                        nepath[t].push_back(i);
                    }
                    path[i][t] = min(path[i][t], path[j][t]);
                }
                for (k = 0; k < nepath[j].ls; k++) {
                    int t = nepath[j].adj[k];
                    if (del[t]) continue;
                    if (path[t][i] == INF) {
                        epath[t].push_back(i);
                        nepath[i].push_back(t);
                    }
                    if (path[t][i] > path[t][j] - path[pre[j]][j]) {
                        path[t][i] = path[t][j] - path[pre[j]][j];
                        fold[i] = j;
                        fpre[i] = j;
                    }
                }
            }
            break;
        }
        if (i > n) {
            for (i = 1; i <= n; i++) {
                if (del[i] || i == root) continue;
                min_cost += path[pre[i]][i];
            }
            break;
        }
    }
    return true;
}

第六章， 几何算法
/*
叉乘
两个点的距离
点到直线距离
返回直线Ax + By + C = 0 的系数
线段
圆
两个圆的公共面积
矩形
根据下标返回多边形的边
两个矩形的公共面积
多边形，逆时针或顺时针给出x，y
多边形顶点
多边形的边
多边形的周长
判断点是否在线段上
判断两条线段是否相交，端点重合算相交
判断两条线段是否平行
判断两条直线是否相交
直线相交的交点
判断是否简单多边形
求多边形面积
判断是否在多边形上
判断是否在多边形内部
点阵的凸包，返回一个多边形
最近点对距离
*/

#include <cmath>
#include <cstdio>
#include <memory>
#include <algorithm>
#include <iostream>
using namespace std;

typedef double TYPE;

#define Abs(x) (((x)>0)?(x):(-(x)))
#define Sgn(x) (((x)<0)?(-1):(1))
#define Max(a,b) (((a)>(b))?(a):(b))
#define Min(a,b) (((a)<(b))?(a):(b))

#define Epsilon 1e-10
#define Infinity 1e+10
#define PI 3.14159265358979323846

TYPE Deg2Rad(TYPE deg) { return (deg * PI) / 180.0;}
TYPE Rad2Deg(TYPE rad) { return (rad * 180.0 / PI);}
TYPE Sin(TYPE deg) {return sin(Deg2Rad(deg));}
TYPE Cos(TYPE deg) { return cos(Deg2Rad(deg));}
TYPE ArcSin(TYPE val) {return Rad2Deg(asin(val));}
TYPE ArcCos(TYPE val) {return Rad2Deg(acos(val));}
TYPE Sqrt(TYPE val) { return sqrt(val);}
struct POINT {
    TYPE x, y, z;
    POINT():x(0), y(0), z(0) {}
    POINT(TYPE _x_, TYPE _y_, TYPE _z_ = 0):
    x(_x_), y(_y_), z(_z_) {}
};
// cross product of (o->a) and (o->b)
// 叉乘
TYPE Cross(const POINT& a, const POINT& b, const POINT& o) {
    return (a.x - o.x) * (b.y - o.y) - (b.x - o.x) * (a.y - o.y);
}

// planar points' distance
// 两个点的距离
TYPE Distance(const POINT& a, const POINT& b) {
    return Sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

struct LINE {
    POINT a;
    POINT b;
    LINE() {}
    LINE(POINT _a_, POINT _b_) : a(_a_), b(_b_) {}
};

// 点到直线距离
double PointToLine(POINT p0, POINT p1, POINT p2, POINT &cp) {
    double d = Distance(p1, p2);
    double s = Cross(p1, p2, p0) / d;
    cp.x = p0.x + s * (p2.y - p1.y) / d;
    cp.y = p0.y - s * (p2.x - p1.x) / d;
    return Abs(s);
}

// 返回 Ax + By + C = 0 的系数
void Coefficient(const LINE& L, TYPE& A, TYPE& B, TYPE& C) {
    A = L.b.y - L.a.y;
    B = L.a.x - L.b.x;
    C = L.b.x * L.a.y - L.a.x * L.b.y;
}

void Coefficient(const POINT& p, const TYPE a, TYPE& A, TYPE& B, TYPE& C) {
    A = Cos(a);
    B = Sin(a);
    C = - (p.y * B + p.x * A);
}

// 线段
struct SEG {
    POINT a;
    POINT b;
    SEG() {}
    SEG(POINT _a_, POINT _b_):a(_a_), b(_b_) {}
};

// 圆
struct CIRCLE {
    TYPE x;
    TYPE y;
    TYPE r;
    CIRCLE() {}
    CIRCLE(TYPE _x_, TYPE _y_, TYPE _r_):x(_x_), y(_y_),r(_r_){}
};

POINT Center (const CIRCLE& circle) {
    return POINT(circle.x, circle.y);
}


TYPE Area(const CIRCLE& circle) {
    return PI * circle.r * circle.r;
}

// 两个圆的公共面积
TYPE CommonArea(const CIRCLE& A, const CIRCLE& B) {
    TYPE area = 0.0;
    const CIRCLE& M = (A.r > B.r) ? A : B;
    const CIRCLE& N = (A.r > B.r) ? B : A;
    TYPE D = Distance(Center(M), Center(N));
    if ((D < M.r + N.r) && (D > M.r - N.r)) {
        TYPE cosM = (M.r * M.r + D * D - N.r * N.r) / (2.0 * M.r * D);
        TYPE cosN = (N.r * N.r + D * D - M.r * M.r) / (2.0 * N.r * D);

        TYPE alpha = 2.0 * ArcCos(cosM);
        TYPE beta = 2.0 * ArcCos(cosN);

        TYPE TM = 0.5 * M.r * M.r * Sin(alpha);
        TYPE TN = 0.5 * N.r * N.r * Sin(beta);

        TYPE FM = (alpha / 360.0) * Area(M);
        TYPE FN = (beta / 360.0) * Area(N);

        area = FM + FN - TM - TN;
    } else if (D <= M.r - N.r) {
        area = Area(N);
    }
    return area;
}

// 矩形
// 矩形的线段
struct RECT {
    POINT a;
    POINT b;
    RECT() {}
    RECT(const POINT& _a_, const POINT & _b_) {a = _a_; b = _b_;}
};

// 根据下标返回多边形的边
SEG Edge(const RECT& rect, int idx) {
    SEG edge;
    while(idx < 0) idx += 4;
    switch(idx % 4) {
        case 0:
            edge.a = rect.a;
            edge.b = POINT(rect.b.x, rect.a.y);
        break;
        case 1:
            edge.a = POINT(rect.b.x, rect.a.y);
            edge.b = rect.b;
        break;
        case 2:
            edge.a = rect.b;
            edge.b = POINT(rect.a.x, rect.b.y);
        break;
        case 3:
            edge.a = POINT(rect.a.x, rect.b.y);
            edge.b = rect.a;
        break;
        default:
        break;
    }
    return edge;
}

TYPE Area(const RECT& rect) {
    return (rect.b.x - rect.a.x) * (rect.b.y - rect.a.y);
}

// 两个矩形的公共面积
TYPE CommonArea(const RECT & A, const RECT & B) {
    TYPE area = 0.0;
    POINT LL(Max(A.a.x, B.a.x), Max(A.a.y, B.a.y));
    POINT UR(Min(A.b.x, B.b.x), Min(A.b.y, B.b.y));

    if((LL.x <= UR.x) && (LL.y <= UR.y)) {
        area = Area(RECT(LL, UR));
    }
    return area;
}

// 多边形，逆时针或顺时针给出x, y
struct POLY {
    int n;
    TYPE* x;
    TYPE* y;
    POLY():n(0), x(NULL), y(NULL) {}
    POLY(int _n_, const TYPE* _x_, const TYPE* _y_) {
        n = _n_;
        x = new TYPE[n + 1];
        memcpy(x, _x_, n * sizeof(TYPE));
        x[n] = _x_[0];

        y = new TYPE[n + 1];
        memcpy(y, _y_, n * sizeof(TYPE));
        y[n] = _y_[0];
    }
};

// 多边形顶点
POINT Vertex(const POLY& poly, int idx) {
    idx %= poly.n;
    return POINT(poly.x[idx], poly.y[idx]);
}

// 多边形的边
SEG Edge(const POLY& poly, int idx) {
    idx %= poly.n;
    return SEG(POINT(poly.x[idx], poly.y[idx]), POINT(poly.x[idx + 1], poly.y[idx + 1]));
}

// 多边形的周长
TYPE Perimeter(const POLY& poly) {
    TYPE p = 0.0;
    for (int i = 0; i < poly.n; i++)
        p = p + Distance(Vertex(poly, i), Vertex(poly, i + 1));
    return p;
}

bool IsEqual(TYPE a, TYPE b) {
    return (Abs(a - b) < Epsilon);
}

bool IsEqual(const POINT& a, const POINT& b) {
    return (IsEqual(a.x, b.x) && IsEqual(a.y, b.y));
}

bool IsEqual(const LINE& A, const LINE& B) {
    TYPE A1, B1, C1;
    TYPE A2, B2, C2;

    Coefficient(A, A1, B1, C1);
    Coefficient(B, A2, B2, C2);

    return IsEqual(A1 * B2 , A2 * B1) &&
            IsEqual(A1 * C2, A2 * C1) &&
            IsEqual(B1 * C2, B2 * C1);
}

// 判断点是否在线段上
bool IsOnSeg(const SEG& seg, const POINT& p) {
    return (IsEqual(p, seg.a) || IsEqual(p, seg.b) ||
            (((p.x - seg.a.x) * (p.x - seg.b.x) < 0 ||
              (p.y - seg.a.y) * (p.y - seg.b.y) < 0 &&
              (IsEqual(Cross(seg.b, p, seg.a), 0))));
}


// 判断两条线段是否相交， 端点重合算相交
bool IsIntersect(const SEG& u, const SEG& v) {
    return (Cross(v.a, u.b, u.a) * Cross(u.b, v.b, u.a) >= 0) &&
            (Cross(u.a, v.b, v.a) * Cross(v.b, u.b, v.a) >= 0) &&
            (Max(u.a.x, u.b.x) >= Min(v.a.x, v.b.x)) &&
            (Max(v.a.x, v.b.x) >= Min(u.a.x, u.b.x)) &&
            (Max(u.a.y, u.b.y) >= Min(v.a.y, v.b.y)) &&
            (Max(v.a.y, v.b.y) >= Min(u.a.y, u.b.y));
}

// 判断两条线段是否平行
bool IsParallel(const LINE& A, const LINE & B) {
    TYPE A1, B1, C1;
    TYPE A2, B2, C2;

    Coefficient(A, A1, B1, C1);
    Coefficient(B, A2, B2, C2);

    return (A1 * B2 == A2 * B1) &&
            ((A1 * C2 != A2 * C1) || (B1 * C2 != B2 * C1))
}

// 判断两条直线是否相交
bool IsIntersect(const LINE& A, const LINE& B) {
    return !IsParallel(A, B);
}

// 直线相交的交点
POINT Intersection(const LINE& A, const LINE& B) {
    TYPE A1, B1, C1;
    TYPE A2, B2, C2;

    Coefficient(A, A1, B1, C1);
    Coefficient(B, A2, B2, C2);

    POINT I(0, 0);
    I.x = -(B2 * C1 - B1 * C2) / (A1 * B2 - A2 * B1);
    I.y = (A2 * C1 - A1 * C2) / (A1 * B2 - A2 * B1);
    return I;
}

bool IsInCircle(const CIRCLE& circle, const RECT& rect) {
    return (circle.x - circle.r >= rect.a.x) &&
            (circle.x + circle.r <= rect.b.x) &&
            (circle.y - circle.r >= rect.a.y) &&
            (circle.y + circle.r <= rect.b.y);
}

// 判断是否简单多边形
bool IsSimple(const POLY& poly) {
    if (poly.n < 3) return false;
    SEG L1, L2;
    for (int i = 0; i < poly.n - 1; i++) {
        L1 = Edge(poly, i);
        for (int j = i + 1; j < poly.n; j++) {
            L2 = Edge(poly, i);
            if (j == i + 1) {
                if (IsOnSeg(L1, L2.b) || IsOnSeg(L2, L1.a)) return false;
            } else if (j == poly.n - i - 1) {
                if (IsOnSeg(L1, L2.a) || IsOnSeg(L2, L1.b)) return false;
            } else {
                if (IsIntersect(L1, L2)) return false;
            }
        }
    }
    return true;
}
// 求多边形面积
TYPE Area(const POLY& poly) {
    if (poly.n < 3) return TYPE(0);
    double s = poly.y[0] * (poly.x[poly.n - 1] - poly.x[1]);
    for (int i = 1; i < poly.n; i++) {
        s += poly.y[i] * (poly.x[i - 1] - poly.x[(i + 1) % poly.n]);
    }
    return s/2;
}

// 判断是否在多边形上
bool IsOnPoly(const POLY & poly, const POINT& p) {
    for (int i = 0; i < poly.n; i++) {
        if (IsOnSeg(Edge(poly, i), p)) return true;
    }
    return false;
}

//判断是否在多边形内部
bool IsInPoly(const POLY& poly, const POINT& p) {
    SEG L(p, POINT(Infinity, p.y));
    int count = 0;
    for (int i = 0; i < poly.n; i++) {
        SEG S = Edge(poly, i);
        if (IsOnSeg(S, p)) return false;
        if (!IsEqual(S.a.y, S.b.y)) {
            POINT& q = (S.a.y > S.b.y) ? (S.a) : (S.b);
            if (IsOnSeg(L, q))  ++count;
            else if (!IsOnSeg(L, S.a) && !IsOnSeg(L, S.b) && IsIntersect(S, L)) ++count;
        }
    }
    return (count % 2 != 0);
}

// 点阵的凸包， 返回一个多边形
POLY ConvexHull(const POINT* set, int n) { // 不适用于点少于3个的情况
    POINT* points = new POINT[n];
    memcpy(points, set, n * sizeof(POINT));

    TYPE* X = new TYPE[n];
    TYPE* Y = new TYPE[n];

    int i, j, k = 0, top = 2;

    for (i = 1; i < n; i++) {
        if ((points[i].y < points[k].y) ||
            ((points[i].y == points[k].y) &&
             (points[i].x < points[k].x))) k = i;
    }

    std::swap(points[0], points[k]);

    for (i = 1; i < n - 1; i++) {
        k = i;
        for (j = i + 1; j < n; j++) {
            if ((Cross(points[j], points[k], points[0]) > 0) ||
                ((Cross(points[j], points[k], points[0]) == 0) &&
                 (Distance(points[0], points[j]) < Distance(points[0], points[k])))) k = j;
        }
        std::swap(points[i], points[k]);
    }

    X[0] = points[0].x; Y[0] = points[0].y;
    X[1] = points[1].x; Y[1] = points[1].y;
    X[2] = points[2].x; Y[2] = points[2].y;

    for (i = 3; i < n; i++) {
        while(Cross(points[i], POINT(X[top], Y[top]),
                    POINT(X[top - 1], Y[top - 1])) >= 0 && top > 0) top--;
        ++top;
        X[top] = points[i].x;
        Y[top] = points[i].y;
    }

    delete [] points;
    POLY poly(++top, X, Y);
    delete [] X;
    delete [] Y;
    return poly;
}

// 最近点对的距离， Written By PrincessSnow
#define MAXN 100000
POINT pt[MAXN];

bool cmp(POINT n1, POINT n2) {
    return (n1.x < n2.x || n1.x == n2.x && n1.y < n2.y);
}

double Get(double dis, int mid, int start, int end) {
    int s = mid, e = mid, i, j;
    double t;

    while(s > start && pt[mid].x - pt[s].x <= dis) s--;
    while(e < end && pt[e].x - pt[mid].x <= dis) e++;
    for (i = s; i <= e; i++)
        for (j = i + 1; j <= e && j <= i + 7; j++) {
            t = Distance(pt[i], pt[j]);
            if (t < dis) dis = t;
        }
    return dis;
}

double ClosestPairDistance(int start, int end) {
    int m = end - start + 1, mid , i;
    double t1, t2, dis = -1, t ;
    if (m <= 3) {
        for (i = start; i < end; i++) {
            t = Distance(pt[i], pt[i + 1]);
            if (t < dis || dis == -1) dis = t;
        }
        t = Distance(pt[start], pt[end]);
        if (t < dis) dis = t;
        return dis;
    }

    if (m % 2 == 0) mid = start + m/2 - 1;
    else mid = start + m/2;
    if (m % 2 == 0) {
        t1 = ClosestPairDistance(start, mid);
        t2 = ClosestPairDistance(mid + 1, end);
    }
    if (t1 < t2) dis = t1;
    else dis = t2;
    dis = Get(dis, mid, start, end);
    return dis;
}

1. 球面上两点最短距离
// 计算圆心角lat 表示纬度， -90 <= w <= 90, lng 表示经度
// 返回两点所在大圆劣弧对应圆心角， 0 <= angle <= pi
double angle(double lng1, double lat1, double lng2, double lat2) {
    double dlng = fabs(lng1 - lng2) * pi / 180;
    while(dlng >= pi + pi) dlng -= pi + pi;
    if (dlng > pi) dlng = pi + pi - dlng;
    lat1 *= pi / 180, lat2 *= pi / 180;
    return acos(cos(lat1) * cos(lat2) * cos(dlng) + sin(lat1) * sin(lat2));
}

// 计算距离，r为球半径
double line_dist(double r, double lng1, double lat1, double lng2, double lat2) {
    double dlng = fabs(lng1 - lng2) * pi / 180;
    while(dlng >= pi + pi) dlng -= pi + pi;
    if (dlng > pi) dlng = pi +pi - dlng;
    lat1 *= pi / 180, lat2 *= pi / 180;
    return r * sqrt(2 - 2 * (cos(lat1) * cos(lat2) * cos(dlng) + sin(lat1) * sin(lat2)));
}

// 计算球面距离， r为球半径
double sphere_dist(double r, double lng1, double lat1, double lng2, double lat2) {
    return r * angle(lng1, lat1, lng2, lat2);
}

2. 三点求圆心坐标
double GetRadiusBy3Points(double x1, double y1,
                          double x2, double y2,
                          double x3, double y3,
                          double &x, double &y) {

    double a11, a12, a21, a22, b1, b2;
    double d, d1, d2;
    a11 = 2 * (x3 - x2);
    a12 = 2 * (y3 - y2);
    a21 = 2 * (x2 - x1);
    a22 = 2 * (y2 - y1);

    b1 = x3 * x3 - x2 * x2 + y3 * y3 - y2 * y2;
    b2 = x2 * x2 - x1 * x1 + y2 * y2 - y1 * y1;

    d = a11 * a22 - a12 * a21;
    d1 = b1 * a22 - a12 * b2;
    d2 = a11 * b2 - b1 * a21;

    x = d1/ d;
    y = d2 / d;
    return (x1 - x) * (x1 - x) + (y1 - y) * (y1 - y);
}

3. 三角形几个重要的点
设三角形的三条边为a，b，c，且不妨设a <= b <= c
三角形的面积可以根据海伦公式：
s = sqrt(p * (p - a) * (p - b) * (p - c)), p = (a + b + c) / 2;

1) 费马点（该点到三角形三个顶点的距离之和最小）
有个有趣的结论： 若三角形的三个内角均小于120，
那么该点连接3个顶点形成的三个角均为120度；
若三角形存在一个内角大于120度，则该顶点就是费马点。
计算公式如下：
若有一个内角大于120度（这里假设为角C）， 则距离为a + b
若三个内角均小于120度， 则距离为
sqrt((a * a + b * b + c * c + 4 * sqrt(3.0) * s) / 2);

2) 内心 -- 角平分线交点
令 x = (a + b - c) / 2, y = (a - b + c) / 2, z = (-a + b + c) / 2,
h = s / p 计算公式为：
sqrt(x * x  + h * h) + sqrt(y * y + h * h) + sqrt(z * z + h * h);

3. 重心 -- 中线交点
计算公式如下：
2.0 / 3 * (sqrt((2 * (a * a + b * b) - c * c) / 4)
           + sqrt((2 * (a * a + c * c) - b * b) / 4)
           + sqrt((2 * (b * b + c * c) - a * a) / 4));
4. 垂心 --- 垂线的交点
计算公式如下：
3 * (c / 2 / sqrt(1 - cosC * cosC))

第七章 专题讨论
1. 树状数组

#include <cstdio>
int data[50001], s[50001], T[50001];
inline int lowbit(int t) {
    return t & (-t);
}
inline int sum(int end) {
    int sum = 0;
    while(end > 0) {
        sum += T[end];
        end -= lowbit(end);
    }
    return sum;
}
inline void plus(int pos, int num, int count) {
    while(pos <= count) {
        T[pos] += num;
        pos += lowbit(pos);
    }
}

int main() {
    char buffer[10];
    int i, j, t, n, a, b;
    scanf("%d", &t);
    for (i = 1; i <= t; i++) {
        scanf("%d", &n);
        T[0] = s[0] = data[0] = 0;
        for (j = 1; j <= n; j++) {
            scanf("%d", &data[j]);
            s[j] = s[j - 1] + data[j];
            T[j] = s[j] - s[j - lowbit(j)];
        }
    }
    return 0;
}
