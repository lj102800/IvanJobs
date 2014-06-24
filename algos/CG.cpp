const double PI = 3.14159265358979323846;
const double eps = 0.0000001;
#define Angle(r) ((r)*180.0/PI)
#define Radian(a) ((a) * PI / 180.0)

struct Point {
    Point(){}
    Point(double _x, double _y):x(_x), y(_y){}
    double x, y;
};
double det(double x1, double y1, double x2, double y2) {
    return x1 * y2 - x2 * y1;
}
double cross(Point& o, Point& a, Point& b) {
    return det(a.x - o.x, a.y - o.y, b.x - o.x, b.y - o.y);
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

