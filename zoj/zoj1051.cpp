#include <cstdio>

char cmap[4] = {'.', '!', 'X', '#'};
int dir[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};

int n = 0;
int D[16]={0};
int culture[20][20]={0};
int k[20][20]={0};

void makek() {
	int i, j, m;
	for (i = 0; i < 20; ++i){
		for (j = 0; j < 20; ++j) {
			// culture[i][j]
			k[i][j] = culture[i][j];
			for (m = 0; m < 4; ++m){
				int newi = i + dir[m][0];
				int newj = j + dir[m][1];
				if (newi >= 0 && newi <20 && newj >= 0 && newj < 20)
					k[i][j] += culture[newi][newj];
			}
		}
	}
}

void makec() {
    int i, j;
    for (i = 0; i < 20; ++i) {
        for (j = 0; j < 20; ++j) {
          culture[i][j] += D[k[i][j]];
          if (culture[i][j] > 3) culture[i][j] = 3;
          else if (culture[i][j] < 0) culture[i][j] = 0;
        }    
    }
}    

void display() {
	int i, j;
	for (i = 0; i < 20; ++i){
		for (j = 0; j < 20; ++j) {
			printf("%c", cmap[culture[i][j]]);
		}
		printf("\n");
	}
}


int main() {
	int cs;
	scanf("%d", &cs);
	int c;
	for (c = 1; c <= cs; ++c) {
		scanf("%d", &n);
		int i;
		for (i = 0; i < 16; ++i)
			scanf("%d", &D[i]);

		int j;
		for (i = 0; i < 20; ++i) {
			for (j = 0; j < 20; ++j) {
				scanf("%d", &culture[i][j]);
			}
		}
		
		for (i = 0; i < n; ++i){
			makek();
			makec();
		}						

		display();
		if (c != cs)
			printf("\n");

	}


	return 0;
}