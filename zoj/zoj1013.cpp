#include <cstdio>
#define MAX 510
int n,
	w1, s1, d1,
	w2, s2, d2,
	w3, s3, d3,
	c1, c2, c3, d4;

int carw[MAX]; // weight limit of every caravan
int cars[MAX]; // size limit of every caravan

struct Item {
	 int x, y, z;
};

#define MAXC 510*510

struct Item Beaker0[MAXC];
struct Item Beaker1[MAXC];


int make(){
	// calculate the fist caravan.
	// max number of equipment 1, caravan 1 can load.
	int n1 = (carw[1]/w1) < (cars[1]/s1) ? (carw[1]/w1) : (cars[1]/s1);

	int n2 = (carw[2]/w2) < (cars[2]/s2) ? (carw[2]/w2) : (cars[2]/s2);

	int i, j, k, m, x, y, p, q, dw, ds, a, b ,c;
	int needle0 = 0;
	int needle1 = 0;

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

				// here, we got a, b, c
				needle0++;
				Beaker0[needle0].x = i;
				Beaker0[needle0].y = j;
				Beaker0[needle0].z = c;				
			}
		}
	}
	int idx = 0; // use this var to record the final res idx.
	for (i = 2; i <= n; ++i){
		// we start calculating prefix i caravans
		// if i % 2 == 0, then we use Beaker0 to get Beaker1
		// if i % 2 == 1, then we use Beaker1 to get Beaker0, and overlap Beaker0's data.
		if (i % 2 == 0){
			needle1 = 0;
			idx = 1;
		}
		else{
			needle0 = 0;
			idx = 0;
		}

		n1 = (carw[i]/w1) < (cars[i]/s1) ? (carw[i]/w1) : (cars[i]/s1);
		n2 = (carw[i]/w2) < (cars[i]/s2) ? (carw[i]/w2) : (cars[i]/s2);
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

					// (a, b, c) -> (j, k , c)
					// check if j, k already in Beaker1 or 0?
					int newa, newb, newc, exi;
					if (i % 2 == 0){
						for (p = 1; p <= needle0; ++p){
							newa = Beaker0[p].x + j;
							newb = Beaker0[p].y + k;
							newc = Beaker0[p].z + c;
							
							exi = -1;
							for (q = 1; q <= needle1; ++q){
								if (Beaker1[q].x == newa && Beaker1[q].y == newb){
									exi = q;
									break;
								}
							}

							if (exi != -1){
								Beaker1[exi].z = (Beaker1[exi].z < newc) ? newc : Beaker1[exi].z;
							} else {
								needle1++;
								Beaker1[needle1].x = newa;
								Beaker1[needle1].y = newb;
								Beaker1[needle1].z = newc;
							}
						}	
					
					} else {
						for (p = 1; p <= needle1; ++p){
							newa = Beaker1[p].x + j;
							newb = Beaker1[p].y + k;
							newc = Beaker1[p].z + c;
							
							exi = -1;
							for (q = 1; q <= needle0; ++q){
								if (Beaker0[q].x == newa && Beaker0[q].y == newb){
									exi = q;
									break;
								}
							}

							if (exi != -1){
								Beaker0[exi].z = (Beaker0[exi].z < newc) ? newc : Beaker0[exi].z;
							} else {
								needle0++;
								Beaker0[needle0].x = newa;
								Beaker0[needle0].y = newb;
								Beaker0[needle0].z = newc;
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
		for (i = 1; i <= needle0; ++i){
			a = Beaker0[i].x;
			b = Beaker0[i].y;
			c = Beaker0[i].z;
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
		for (i = 1; i <= needle1; ++i){
			a = Beaker1[i].x;
			b = Beaker1[i].y;
			c = Beaker1[i].z;
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