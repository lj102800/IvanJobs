typedef complex<double> Point;
#define x real()
#define y imag()

double cross(const Point& a, const Point& b) {
	return a.x * b.y - a.y * b.x;
}

bool cmp(const Point& p1, const Point& p2) {
	return arg(p1) < arg(p2);
}

bool cmp(const Point& p1, const Point& p2) {
	return atan2(p1,y, p1.x) < atan2(p2.y, p2.x);
}

bool cmp(const Point& p1, const Point& p2) {
	if (p1.y == 0 && p2.y == 0 && p1.x * p2.x <= 0) return p1.x > p2.x;
	if (p1.y == 0 && p1.x >= 0 && p2.y != 0) return true;
	if (p2.y == 0 && p2.x >= 0 && p1.y != 0) return false;
	if (p1.y * p2.y < 0) return p1.y > p2.y;
	double c = cross(p1, p2);
	return c > 0 || c == 0 && fabs(p1.x) < fabs(p2.x);
}

bool cmp(const Point& p1, const Point& p2) {
	if (p1.y > 0 && p2.y > 0) 
		return p2.x * p1.y < p2.y * p1.x;
	else if (p1.y < 0 && p2.y < 0)
		return p2.x * p1.y < p2.y * p1.x;
	else if (p1.y == 0)
		if (p1.x > 0) return true;
		else return p2.y < 0;
	else if (p2.y == 0)
		if (p2.x > 0) return false;
		else return p1.y > 0;
	else return p1.y > 0;
}

void SortPointsByPolarAngle() {
	Point p[100];
	sort(p, p + 100, cmp);
}
