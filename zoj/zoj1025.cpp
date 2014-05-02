#include <cstdio>
#include <cstring>
#include <cstdlib>

#define MAX 5100

struct stick {
	int wei,len;
};

stick sticks[MAX];
int grp[MAX];

int cmp(const void *a, const void *b){
	stick *pre = (stick*)a;
	stick *post = (stick*)b;
	if (pre->len == post->len)
		return pre->wei - post->wei;
	return pre->len - post->len;
}

int main(){
	int T, i, n, j, k;
	scanf("%d", &T);
	for (i = 1; i <= T; ++i){
		scanf("%d", &n);
		for (j = 0; j < n; ++j){
			scanf("%d%d", &(sticks[j].len), &(sticks[j].wei));
		}
			
		qsort(sticks, n, sizeof(stick), cmp);
		memset(grp, 0, sizeof(grp));
		grp[0] = 1;
	
		for (j = 1; j < n; ++j){
			int tmp = 0;
			for (k = 0; k < j; ++k){
				if (sticks[j].wei < sticks[k].wei && tmp < grp[k]) tmp = grp[k]; 
			}
			grp[j] = tmp + 1;
		}

		int res = -1;
		for (j = 0; j < n; ++j){
			if (grp[j] > res)
				res = grp[j];
		}
		printf("%d\n", res);
	}
	return 0;
}