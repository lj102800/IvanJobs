#include <cmath>
const double PI = acos(-1);
const double eps = 1e-6;
#define Angle(r) ((r)*180.0/PI)
#define Radian(a) ((a) * PI / 180.0)
#define Mn(a, b) ((a) < (b) ? (a) : (b))
#define Mx(a, b) ((a) > (b) ? (a) : (b))
struct Point {
    Point(){}
    Point(double _x, double _y):x(_x), y(_y){}
    double x, y;
};

double cross(Point& o, Point& a, Point& b) {
    double x1 = a.x - o.x; double y1 = a.y - o.y;
    double x2 = b.x - o.x; double y2 = b.y - o.y;
    return x1 * y2 - x2 * y1;
}
double dot(Point& o, Point& a, Point& b) {
    double x1 = a.x - o.x; double y1 = a.y - o.y;
    double x2 = b.x - o.x; double y2 = b.y - o.y;
    return x1 * x2 + y1 * y2;
}
double dis(Point& a, Point& b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return sqrt(dx * dx + dy * dy);
}
int dblcmp(double v1, double v2) {
    double delta = v1 - v2;
    if (fabs(delta) < eps) return 0;
    return (delta > 0) ? 1 : -1;
}
bool inter(Point& a, Point& b, Point& c, Point& d) {
    if (Mn(a.x, b.x) > Mx(c.x, d.x) ||
	Mn(a.y, b.y) > Mx(c.y, d.y) ||
	Mn(c.x, d.x) > Mx(a.x, b.x) ||
	Mn(c.y, d.y) > Mx(a.y, b.y)) return false;
    double h, i, j, k;
    h = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    i = (b.x - a.x) * (d.y - a.y) - (b.y - a.y) * (d.x - a.x);
    j = (d.x - c.x) * (a.y - c.y) - (d.y - c.y) * (a.x - c.x);
    k = (d.x - c.x) * (b.y - c.y) - (d.y - c.y) * (b.x - c.x);
    return dblcmp(h * i, 0.0) <= 0 && dblcmp(j * k, 0.0) <= 0;
}


