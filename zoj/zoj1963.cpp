#include <cstdio>
#include <cstring>

int seconds(int h, int m, int s){
	return s + m*60 + 60*60*h;
}

char timestr[25];

int main(){
	// n, the number of sections;  1<=n<=20
	// d, total distance of the relay; 0.0 < d < 200.0
	int n, t;
	float d;
	scanf("%d", &n);
	scanf("%f", &d);
	while(scanf("%d", &t) != EOF){
		int i;
		int h, m, s;
		int resm, ress;
		bool isdis = false;
		float totals = 0.0;
		for (i = 0; i < n; ++i){
			scanf("%s", timestr);	
			if (strcmp(timestr, "-:--:--") == 0)
				isdis = true;
			else{
				sscanf(timestr, "%d:%d:%d", &h, &m, &s);
				totals += seconds(h, m, s);
			}
		}
		float te1 = totals / d;
		int te2 = (int)(te1 + 0.5);
		ress = te2 % 60;
		resm = te2 / 60;

		if(!isdis){
			printf("%3d: %d:%02d min/km\n", t, resm, ress);
		} else {
			printf("%3d: -\n", t);
		}
	}
	return 0;
}