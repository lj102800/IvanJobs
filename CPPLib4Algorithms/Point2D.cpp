class Point2D {
	private:
		double x;
		double y;
	public:
		Point2D(double _x, double _y) : x(_x), y(_y) {}
		// whether a b c form a counter-clockwise 
		static int ccw(Point2D& a, Point2D& b, Point2D& c) {
			double area2 = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
			if (area2 < 0) return -1; // clockwise
			else if (area2 > 0) return 1; // counter-clockwise
			else return 0; // colinear
		}
};
