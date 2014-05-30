const double eps = 1e-8;

int dcmp(double d) {
	return d < -eps ? -1 : d > eps;
}

struct Point {
    double x, y;

    Point() {}
    Point(double x, double y)
    : x(x), y(y) {}

    Point operator + (const Point& a) const {
        return Point(x + a.x, y + a.y);
    }
    Point operator - (const Point& a) const {
        return Point(x - a.x, y - a.y);
    }
    Point operator * (double k) const {
        return Point(x * k, y * k);
    }
    Point operator / (double k) const {
        return Point(x / k, y / k);
    }

    bool operator < (const Point &a) {
        return x < a.x || (x == a.x && y < a.y);
    }
    bool operator == (const Point &a) {
        return dcmp(x - a.x) == 0 && dcmp(y - a.y) == 0;
    }
};
double dot(const Point &a, const Point &b) {
    return a.x * b.x + a.y * b.y;
}
double cross(const Point &a, const Point &b) {
    return a.x * b.y - a.y * b.x;
}
double length(const Point &a) {
    return sqrt(dot(a, a));
}
