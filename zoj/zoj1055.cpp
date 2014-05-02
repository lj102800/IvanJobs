#include <cstdio>
#include <cstring>
#include <deque>

using namespace std;

struct LOCATION {
	int x, y;
};

#define MATRIXLEN 30
#define LETTERLEN 26
#define INFINITE 2<<30
char Mat[MATRIXLEN][MATRIXLEN];
int D[MATRIXLEN][MATRIXLEN];
char P[MATRIXLEN][MATRIXLEN]; // 'U'->up, 'D'->down, 'L'->left, 'R'->right, '\0' -> NIL
char Color[MATRIXLEN][MATRIXLEN]; // 'W' -> 0, 'B'->1, 'G'->2
double Load[MATRIXLEN][MATRIXLEN];
int COLN, ROWN;
//here, we do a map on (x, y) -> x.y
LOCATION Door[LETTERLEN];
char Dir[][2]={{0, -1}, {0, 1}, {-1, 0}, {1, 0}};
char DirL[] = {'L', 'R', 'U', 'D'};

deque<LOCATION> Q;

void rec(int x, int y, char d, double load){
	if (Mat[x][y] == d)
		return ;
	Load[x][y] += load;
	int du = D[x][y];
	int i, nx, ny, count;
	count = 0;
	for (i = 0; i < 4; ++i){
		nx = x + Dir[i][0];
		ny = y + Dir[i][1];
		if (nx >= 0 && nx < ROWN && ny >= 0 && ny < COLN){
			if (Mat[nx][ny] == '.' && D[nx][ny] == du - 1){
				count++;
			}	
		}
	}
	if (count == 0)
		return ;
	else {
		load = load / count;
		for (i = 0; i < 4; ++i){
			nx = x + Dir[i][0];
			ny = y + Dir[i][1];
			if (nx >= 0 && nx < ROWN && ny >= 0 && ny < COLN){
				if (Mat[nx][ny] == '.' && D[nx][ny] == du - 1){
					rec(nx, ny, d, load);
				}	
			}
		}
		load = load * count;
	}
}

void make(char s, char d, double load){
	// base on s, do BFS, get D, P, Color done.
	memset(Color, 0, sizeof(Color));
	memset(P, 0, sizeof(P));
	int i, j;
	for (i = 0; i < MATRIXLEN; ++i){
		for (j = 0; j < MATRIXLEN; ++j){
			D[i][j] = INFINITE;
		}
	}

	int sx , sy;
	int nx, ny;
	int vx, vy;
	//mmap(Door[d - 'A'], vx, vy);
	//mmap(Door[s-'A'], sx, sy);
	vx = Door[d - 'A'].x;
	vy = Door[d - 'A'].y;
	sx = Door[s - 'A'].x;
	sy = Door[s - 'A'].y;
	// here, we get s's x and y
	Color[sx][sy] = 2;
	D[sx][sy] = 0;
	P[sx][sy] = 0;
	
	Q.clear();
	Q.push_back(Door[s - 'A']);
	while(!Q.empty()){
		LOCATION u = Q.front();
		Q.pop_front();
		int ux, uy;
		ux = u.x;
		uy = u.y;	
		for (i = 0; i < 4; ++i){
			nx = ux + Dir[i][0];
			ny = uy + Dir[i][1];
			if (nx >= 0 && nx < ROWN && ny >= 0 && ny < COLN){
				if (Color[nx][ny] == 0 && Mat[nx][ny] != 'X'){
					Color[nx][ny] = 2;
					D[nx][ny] = D[ux][uy] + 1;
				
					switch(DirL[i]){
						case 'L':
							P[nx][ny] = 'R';
							break;
						case 'R':
							P[nx][ny] = 'L';
							break;
						case 'U':
							P[nx][ny] = 'D';
							break;
						case 'D':
							P[nx][ny] = 'U';
							break;
					}
					LOCATION tmp;
					tmp.x = nx;
					tmp.y = ny;
					Q.push_back(tmp);
				}
			}
		}
		Color[ux][uy] = 1;
	}
	
	
	//after BFS
	int count = 0;
	for (i = 0; i < 4; ++i){
		nx = vx + Dir[i][0];
		ny = vy + Dir[i][1];
		if (nx >= 0 && nx < ROWN && ny >= 0 && ny < COLN){
			if (Mat[nx][ny] == '.' && D[nx][ny] == D[vx][vy] - 1){
				count++;
			}
		}
	}
	if (count == 0)
		return ;
	for (i = 0; i < 4; ++i){
		nx = vx + Dir[i][0];
		ny = vy + Dir[i][1];
		if (nx >= 0 && nx < ROWN && ny >= 0 && ny < COLN){
			if (Mat[nx][ny] == '.' && D[nx][ny] == D[vx][vy] - 1){
				rec(nx, ny, d, load / count);
			}
		}
	}
}

int main(){
	
	scanf("%d%d", &COLN, &ROWN);
	int i, j, len;
	for (i = 0; i < ROWN; ++i){
		scanf("%s", Mat[i]);
		len = strlen(Mat[i]);
		for (j = 0; j < len; ++j){
			char tmpc = Mat[i][j];
			if (tmpc != 'X' && tmpc != '.'){
				Door[tmpc - 'A'].x = i;
				Door[tmpc - 'A'].y = j;
			}
		}
	}
	char SD[10];
	int load;
	while(scanf("%s", SD) && !(SD[0] == 'X' && SD[1] == 'X')){
		scanf("%d", &load);
		// here, we got SD, and load
		make(SD[0], SD[1], double(load));
	}

	for (i = 0; i < ROWN; ++i){
		for (j = 0; j < COLN; ++j){
			if (j != 0)
				printf(" ");
			printf("%6.2lf", Load[i][j]);
		}
		printf("\n");
	}
	//system("PAUSE");
	return 0;
}