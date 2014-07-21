/**
 Graham Scan for HDU 1392
 */
#include <cstdio>
#include <algorithm>
#include <stack>
#include <cmath>
using namespace std;

#define REP(i, n) for(int _n = n, i = 0; i < _n; i++)
#define FOR(i, a, b) for(int i = (a), _b = (b); i <= _b; i++)
#define Max(a, b) ((a) > (b) ? (a) : (b))
#define Min(a, b) ((a) < (b) ? (a) : (b))
#define L(fmt, ...) do {if(false) printf(fmt"\n", ##__VA_ARGS__);} while(false)

struct Point{
	Point():x(0.0), y(0.0) {}
	Point(double _x, double _y) : x(_x), y(_y) {}
	double x,y;
};
// cross product of two vectors(represented by struct Point)
double cross(const Point& a, const Point& b) {
	return a.x * b.y - a.y * b.x;
}
// polar angle sort, base unit vector is (1, 0) ? increasing order?
bool cmp(const Point& p1, const Point& p2) {
	if (p1.y == 0 && p2.y == 0 && p1.x * p2.x <= 0) return p1.x > p2.x;
	if (p1.y == 0 && p1.x >= 0 && p2.y != 0) return true;
	if (p2.y == 0 && p2.x >= 0 && p1.y != 0) return false;
	if (p1.y * p2.y < 0) return p1.y > p2.y;
	double c = cross(p1, p2);
	return c > 0 || c == 0 && fabs(p1.x) < fabs(p2.x);
}
// use directional area to indicate whether points a, b, c are counter-clockwise
int ccw(Point& a, Point& b, Point& c) {
	double area2 = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
	if (area2 < 0) return -1; // clockwise
	else if (area2 > 0) return 1; // counter-clockwise
	else return 0; // colinear
}
bool comp(const Point& p1, const Point& p2) {
	if (p1.y < p2.y) return true;
	else if (p1.y == p2.y) return p1.x < p2.x;
	return false;
}
double dis(Point& a, Point& b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return sqrt(dx * dx + dy * dy);
}
#define MAXN 110
Point T[MAXN];
int main() {
	int n;
	while(scanf("%d", &n) && n) {
		double res;
		REP(i, n) scanf("%lf%lf", &(T[i].x), &(T[i].y));
		if (n == 1) {
			printf("0.00\n");
			continue;
		}
		if (n == 2) {
			res = dis(T[0], T[1]);
			printf("%.2lf\n", res);
			continue;
		}
		sort(T, T + n, comp);
		FOR(i, 1, n - 1) {
			T[i].x -= T[0].x;
			T[i].y -= T[0].y;
		}
		sort(T + 1, T + n, cmp);
		FOR(i, 1, n - 1) {
			T[i].x += T[0].x;
			T[i].y += T[0].y;
		}
		stack<Point> sp;
		sp.push(T[0]); sp.push(T[1]);
		T[n].x = T[0].x; T[n].y = T[0].y;
		FOR(i, 2, n) {
			sp.push(T[i]);
			while(sp.size() > 2) {
				Point p3 = sp.top(); sp.pop();
				Point p2 = sp.top(); sp.pop();
				Point p1 = sp.top(); sp.pop();
				L("(%lf, %lf) (%lf, %lf) (%lf, %lf)", p1.x, p1.y, p2.x, p2.y, p3.x, p3.y);
				if (ccw(p1, p2, p3) == -1) {
					L("false");
					sp.push(p1); sp.push(p3);
				} else{
					L("true");
					sp.push(p1); sp.push(p2); sp.push(p3);
				 	break;
				}
			}
		}
		int num = sp.size();
		L("num:%d", num);
		REP(i, num) {
			T[i] = sp.top(); sp.pop();
		}	
		res = 0.0;
		FOR(i, 1, num - 1) {
			res += dis(T[i - 1], T[i]);
		}
		printf("%.2lf\n", res);
	}
	return 0;
}

