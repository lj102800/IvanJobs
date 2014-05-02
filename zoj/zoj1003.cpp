#include <cstdio>
#include <cstring>


int record[110]; // use index 1 - 100 to record, 0 : you can use it, 1 : you can't use it
bool isavalid = false;
bool isbvalid = false;
int a, b;


void judge(int n){
	if (n == 1){
		isbvalid = true;
		return ;
	}

	int i;
	for(i = 1; i <= 100; ++i){
		if (record[i] == 0 && n % i == 0){
			n /= i;
			record[i] = 1;
			judge(n);
			n *= i;
			record[i] = 0;
		}
	}
}

void make(int m){
	if (m == 1){
		isavalid = true;
		// we find one case of m here, judge n now
		judge(b);
		return ;
	}

	int i;
	for (i = 1; i <= 100; ++i){
		if (record[i] == 0 && m % i == 0){
			record[i] = 1;			
			m /= i;
			make(m);
			m *= i;
			record[i] = 0;
		}
	}
}

int main(){
	int m, n;
	while (scanf("%d%d", &m, &n) != EOF){
		a = m < n ? m : n;
		b = m < n ? n : m;
		// a < b
		memset(record, 0, sizeof record);
		isavalid = false;
		isbvalid = false;
		int ra = a;
		int rb = b;

		make(a);

		if (isavalid && !isbvalid){
			printf("%d\n", ra);
		} else {
			printf("%d\n", rb);
		}
	}

	return 0;
}