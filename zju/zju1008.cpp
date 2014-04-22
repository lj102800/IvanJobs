#include <cstdio>
#include <cstring>

struct tri {
	int up;
	int right;
	int bottom;
	int left;
};

tri a[100];
tri b[10][10];
int count[100];

int n;
int cur;
bool possible = false;
bool inrange(int i, int j){
	if (i>=0 && i<=n-1 && j>=0 && j<=n-1)
		return true;
	else
		return false;
}
void make(int i, int j){
	if (i == n && j == 0){
		// we find one
		possible = true;
		return ;
	}

	int m;
	for (m = 0; m < cur; ++m){
		if(count[m] != 0){
			// b[i-1][j] b[i][j-1] vs a[m]
			bool isok = true;
			if(inrange(i-1, j)){
				if(b[i-1][j].bottom != a[m].up)
					isok = false;
			}
			if(inrange(i, j-1)){
				if(b[i][j-1].right != a[m].left)
					isok = false;
			}
			if(isok){
				b[i][j].up = a[m].up;
				b[i][j].bottom = a[m].bottom;
				b[i][j].left = a[m].left;
				b[i][j].right = a[m].right;
				count[m]--;

				if(j == n-1 && !possible)
					make(i+1, 0);
				if(j != n-1 && !possible)
					make(i, j+1);

				// backtrace
				count[m]++;
			}
		}
	}
}

int main(){
	int num = 0;
	
	while (scanf("%d", &n) && n!=0){
		int i, j;
		cur = 0;
		memset(count, 0, sizeof count);
		for (i = 0; i < n*n; ++i){
			int u, r, b, lef;
			scanf("%d", &u);
			scanf("%d", &r);
			scanf("%d", &b);
			scanf("%d", &lef);
			bool find = false;
			for (j = 0; j <cur; ++j){
				if(a[j].up == u &&
				   a[j].bottom == b &&
				   a[j].left == lef &&
				   a[j].right == r){
					count[j]++;
					find = true;
				}
			}
			if(!find){
				// if not find
				a[cur].up = u;
				a[cur].right = r;
				a[cur].bottom = b;
				a[cur].left = lef;
				count[cur]++;
				cur++;
			}

		}
		possible = false;
		make(0, 0);
		num++;
		if(num != 1)
				printf("\n");
		if(possible){
			
			printf("Game %d: Possible\n", num);
			
		}
		else{
			printf("Game %d: Impossible\n", num);
		}
	}
	return 0;
}
