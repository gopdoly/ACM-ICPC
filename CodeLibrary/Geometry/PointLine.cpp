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
    void read() {
        scanf("%lf%lf", &x, &y);
    }

    bool operator < (const Point &a) const {
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

double angle(Point a, Point b) {
    return acos(dot(a, b) / length(a) / length(b));
}
Point rotate(Point a, double r) {
    return Point(a.x * cos(r) - a.y * sin(r), a.x * sin(r) + a.y * cos(r));
}
Point normal(Point a) {//left
    return Point(-a.y, a.x);
}
double distToLine(Point p, Point a, Point b) {
    Point v1 = b - a, v2 = p - a;
    return fabs(cross(v1, v2)) / length(v1);
}
Point pointToLine(Point p, Point a, Point b) {
    Point v = b - a;
    return a + v * (dot(v, p - a) / dot(v, v));
}
double distToSeg(Point p, Point a, Point b) {
    if(a == b) return length(p - a);
    Point v1 = b - a, v2 = p - a, v3 = p - b;
    if(dcmp(dot(v1, v2)) < 0) return length(v2);
    else if(dcmp(dot(v1, v3)) > 0) return length(v3);
    else return fabs(cross(v1, v2)) / length(v1);
}
Point pointToSeg(Point p, Point a, Point b) {
    if(a == b) return a;
    Point v1 = b - a, v2 = p - a, v3 = p - b;
    if(dcmp(dot(v1, v2)) < 0) return a;
    else if(dcmp(dot(v1, v3)) > 0) return b;
    else return a + v1 * (dot(v1, p - a) / dot(v1, v1));
}
Point symPoint(Point p, Point a, Point b) {
    return pointToLine(p, a, b) * 2 - p;
}
Point crossPoint(Point a1, Point b1, Point a2, Point b2) {
    double s1 = cross((b2 - a1), (a2 - a1));
    double s2 = cross((b2 - b1), (a2 - b1));
    return (a1 * s2 - b1 * s1) / (s2 - s1);
}
int segCrossSeg(Point a1, Point b1, Point a2, Point b2) {
    int d1 = dcmp(cross(b1 - a2, a1 - a2));
    int d2 = dcmp(cross(b1 - b2, a1 - b2));
    int d3 = dcmp(cross(b2 - a1, a2 - a1));
    int d4 = dcmp(cross(b2 - b1, a2 - b1));
    if ((d1 ^ d2) == -2 && (d3 ^ d4) == -2) return 2;
    return ((d1 == 0 && dcmp(dot(b1 - a2, a1 - a2)) <= 0)
         || (d2 == 0 && dcmp(dot(b1 - b2, a1 - b2)) <= 0)
         || (d3 == 0 && dcmp(dot(b2 - a1, a2 - a1)) <= 0)
         || (d4 == 0 && dcmp(dot(b2 - b1, a2 - b1)) <= 0));
}
int lineCrossSeg(Point a1, Point b1, Point a2, Point b2) {
    int d1 = dcmp(cross(b1 - a2, a1 - a2));
    int d2 = dcmp(cross(b1 - b2, a1 - b2));
    if((d1 ^ d2) == -2) return 2;
    return (d1 == 0 || d2 == 0);
}
bool onSegment(Point p, Point a, Point b) {
    return dcmp(cross(a - p, b - p)) == 0 && dcmp(dot(a - p, b - p) <= 0);
}
