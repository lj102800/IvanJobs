#include <cstdio>
#include <algorithm>

int m, n;
bool make(int m, int n, bool fh){
	int mm = std::max(m, n);
	int nn = std::min(m, n);

	if (mm == nn || mm / 2 >= nn)
		return fh;

	return make(nn, mm - nn, !fh);
}

int main(){
	while(scanf("%d%d", &m, &n) && !(m == 0 && n == 0)){
		bool res = make(m, n, false);
		if (res)
			printf("Ollie wins\n");
		else 
			printf("Stan wins\n");
	}
	return 0;
}
