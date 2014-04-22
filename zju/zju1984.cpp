#include <cstdio>

char mat[5100];
int n;

bool isfind = false;

bool check(int cur){
	int len = cur + 1;
	int i = cur;
	int j;
	while(cur - i + 1 <= len / 2){
		// i..cur
		bool isequal = true;
		for (j = i; j <= cur; ++j){
			if (mat[j] != mat[i - 1 - cur + j]){
				isequal = false;// we find a not equal, 
				break;
			}
		}
		if (isequal) // all char is equal
			return false;
		i--;
	}
	return true;
}

void make(int i){
	if (i == n || isfind){
		isfind = true;
		return ;
	}
	

	mat[i] = 'N';
	if(check(i)){
		make(i + 1);
	}
	
	if (isfind)
		return ;
	mat[i] = 'P';
	if(check(i)){
		make(i + 1);
	}

	if (isfind)
		return ;
	mat[i] = 'O';
	if(check(i)){
		make(i + 1);
	}
}


int main() {
	
	while (scanf("%d", &n) && n != 0){
		mat[0] = 'N';
		isfind = false;

		make(1);
		
		if (isfind){
			int i;
			for (i = 0; i < n; ++i){
				printf("%c", mat[i]);
			}
			printf("\n");
		}else{
			printf("\n");
		}
	}
	return 0;
}
