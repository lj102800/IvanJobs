#include <cstdio>

int board[5][5];
int count = 0;
int n = 0;
int max = -1;

void put(int ith){
	if (ith == n*n){
		max = (max > count ? max : count);
		return ;
	}

	int i, j;
	for (i = ith; i < n*n; ++i){
		int r = i / n ;
		int c = i % n ;
		if (board[r][c] == -1){
			board[r][c] = i;
			count++;
			// mark
			// up
		
			for (j = r-1; j >= 0; j--){
			
				if(board[j][c] == -2)
					break;
				else if (board[j][c] == -1)
					board[j][c] = i;
			}
			// down
			for (j = r + 1; j < n; ++j){
				if(board[j][c] == -2)
					break;
				else if (board[j][c] == -1)
					board[j][c] = i;
			}

			// left
			for (j = c-1; j >= 0; j--){
				if (board[r][j] == -2)
					break;
				else if (board[r][j] == -1)
					board[r][j] = i;
			}

			// right
			for (j = c + 1; j < n; ++j){
				if (board[r][j] == -2)
					break;
				else if (board[r][j] == -1)
					board[r][j] = i;
			}

			put(i + 1);
			count--;

			board[r][c] = -1;
			// up
			for (j = r-1; j >= 0; j--){
				if(board[j][c] == -2)
					break;
				else if (board[j][c] == i)
					board[j][c] = -1;
			}
			// down
			for (j = r + 1; j < n; ++j){
				if(board[j][c] == -2)
					break;
				else if (board[j][c] == i)
					board[j][c] = -1;
			}

			// left
			for (j = c-1; j >= 0; j--){
				if (board[r][j] == -2)
					break;
				else if (board[r][j] == i)
					board[r][j] = -1;
			}

			// right
			for (j = c + 1; j < n; ++j){
				if (board[r][j] == -2)
					break;
				else if (board[r][j] == i)
					board[r][j] = -1;
			}

		} else {
			put(i + 1);
		}
	}
}
int main(){
	while(scanf("%d", &n) && n != 0) {
		getchar();
		int i, j;
		for (i = 0; i < n; ++i){
			for (j = 0; j < n; ++j){
				char c;
				scanf("%c", &c);
				if(c == 'X')
					board[i][j] = -2;
				else
					board[i][j] = -1;
			}
			getchar();
		}

		max = -1;
		count = 0;

		put(0);
		printf("%d\n", max);
	}
	return 0;
}
