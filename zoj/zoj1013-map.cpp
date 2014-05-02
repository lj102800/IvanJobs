#include <cstdio>
#include <map>
using namespace std;

#define MAX 510
int n,
	w1, s1, d1,
	w2, s2, d2,
	w3, s3, d3,
	c1, c2, c3, d4;

int carw[MAX]; // weight limit of every caravan
int cars[MAX]; // size limit of every caravan



#define MAXC 510*510

map<float, int> Beaker0;
map<float, int> Beaker1;


int make(){
	// calculate the fist caravan.
	// max number of equipment 1, caravan 1 can load.
	int n1 = (carw[1]/w1) < (cars[1]/s1) ? (carw[1]/w1) : (cars[1]/s1);

	int n2 = (carw[1]/w2) < (cars[1]/s2) ? (carw[1]/w2) : (cars[1]/s2);

	int i, j, k, m, x, y, p, q, dw, ds, a, b ,c, tmp1, tmp2;
	map<float, int>::iterator it1, it2;
	float first;
	int second;

	Beaker0.clear();
	Beaker1.clear();
	for (i = 0; i <= n1; ++i){
		if (carw[1] - i * w1 < 0 || cars[1] - i * s1 < 0)
			break;
		for (j = 0; j <= n2; ++j){
			dw = carw[1] - (i * w1 + j * w2);
			ds = cars[1] - (i * s1 + j * s2);

			if (dw < 0 || ds < 0)
				break;
			else {
				// calc number of equip 3
				c = 0;
				while (dw >= 0 && ds >= 0){
					c++;
					dw -= w3;
					ds -= s3;
				}
				c -= 1;

				// here, we got i, j, c
				Beaker0[i + float(j)/1000]=c;				
			}
		}
	}
	int idx = 0; // use this var to record the final res idx.
	for (i = 2; i <= n; ++i){
		// we start calculating prefix i caravans
		// if i % 2 == 0, then we use Beaker0 to get Beaker1
		// if i % 2 == 1, then we use Beaker1 to get Beaker0, and overlap Beaker0's data.
		if (i % 2 == 0){
			Beaker1.clear();
			idx = 1;
		}
		else{
			Beaker0.clear();
			idx = 0;
		}
		tmp1 = carw[i]/w1;
		tmp2 = cars[i]/s1;
		n1 = tmp1 < tmp2 ? tmp1 : tmp2;
		tmp1 = carw[i]/w2;
		tmp2 = cars[i]/s2;
		n2 = tmp1 < tmp2 ? tmp1 : tmp2;

		for (j = 0; j <= n1; ++j){
			if (carw[i] - j * w1 < 0 || cars[i] - j * s1 < 0)
				break;
			for (k = 0; k <= n2; ++k){
				dw = carw[i] - (j * w1 + k * w2);
				ds = cars[i] - (j * s1 + k * s2);

				if (dw < 0 || ds < 0){
					break;
				} else {
					// calc number of equip 3
					c = 0;
					while (dw >= 0 && ds >= 0){
						c++;
						dw -= w3;
						ds -= s3;
					}
					c -= 1;

					// (j, k , c)
					// check if j, k already in Beaker1 or 0?
					
					if (i % 2 == 0){
						for (it1 = Beaker0.begin(); it1 != Beaker0.end(); ++it1){
							first = it1->first;
							second = it1->second;
							first += j;
							first += k/1000.0;
							second += c;
							it2 = Beaker1.find(first);
							if (it2 != Beaker1.end()){
								it2->second = it2->second < second ? second : it2->second;
							} else {
								Beaker1[first] = second;
							}
						}					
					} else {
						for (it1 = Beaker1.begin(); it1 != Beaker1.end(); ++it1){
							first = it1->first;
							second = it1->second;
							first += j;
							first += k/1000.0;
							second += c;
							it2 = Beaker0.find(first);
							if (it2 != Beaker0.end()){
								it2->second = it2->second < second ? second : it2->second;
							} else {
								Beaker0[first] = second;
							}
						}					
					}
				}
			}
		}
	}
	
	int res = -1;
	int min, defend;
	if (idx == 0){
		for (it1 = Beaker0.begin(); it1 != Beaker0.end(); ++it1){
			a = int(it1->first);
			b = int((it1->first - a)*1000);
			c = it1->second;
			min = 1000;
			min = a < min ? a : min;
			min = b < min ? b : min;
			min = c < min ? c : min;

			a -= min;
			b -= min;
			c -= min;
				
			defend = 0;
			defend += min * d4;
			defend += a * d1;
			defend += b * d2;
			defend += c * d3;

			res = defend > res ? defend : res;

		}
	} else { // idx == 1
		for (it1 = Beaker1.begin(); it1 != Beaker1.end(); ++it1){
			a = int(it1->first);
			b = int((it1->first - a) * 1000);
			c = it1->second;
			min = 1000;
			min = a < min ? a : min;
			min = b < min ? b : min;
			min = c < min ? c : min;

			a -= min;
			b -= min;
			c -= min;
				
			defend = 0;
			defend += min * d4;
			defend += a * d1;
			defend += b * d2;
			defend += c * d3;

			res = defend > res ? defend : res;

		}
	}
	
	return res;
}

int main(){
	int sh = 0;
	while (scanf("%d", &n) && n != 0){
		sh++;
		if (sh != 1)
			printf("\n");
		scanf("%d%d%d", &w1, &s1, &d1);
		scanf("%d%d%d", &w2, &s2, &d2);
		scanf("%d%d%d", &w3, &s3, &d3);
		scanf("%d%d%d%d", &c1, &c2, &c3, &d4);

		int i;
		for (i = 1; i <= n; ++i){
			scanf("%d%d", &carw[i], &cars[i]);
		}

		int res = make();
		printf("Case %d:%d\n", sh, res);
	}
	return 0;
}