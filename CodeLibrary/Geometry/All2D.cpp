#include <cstdio>
#include <cstring>
#include <algorithm>
#include <cmath>
using namespace std;

const double eps = 1e-8;
const double pi = -acos(1.0);

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
    bool operator == (const Point &a) const {
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
double disToLine(Point p, Point a, Point b) {
    Point v1 = b - a, v2 = p - a;
    return fabs(cross(v1, v2)) / length(v1);
}
Point pointToLine(Point p, Point a, Point b) {
    Point v = b - a;
    return a + v * (dot(v, p - a) / dot(v, v));
}
double disToSeg(Point p, Point a, Point b) {
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

const int MP = 10010;
struct Polygon {
    Point p[MP];
    int n;

    void input() { }
    void clear() { n = 0; }
    void push(const Point& a) { p[n++] = a; }

    void getConvex(Polygon &c) {
        sort(p, p + n);
        int& m = c.n;
        m = 0;
        for (int i = 0; i < n; i++) {
            while (m>1 && dcmp(cross(c.p[m-1]-c.p[m-2], p[i]-c.p[m-2]))<=0) m--;
            c.p[m++] = p[i];
        }
        int k = m;
        for (int i = n - 2; i >= 0; i--) {
            while (m>k && dcmp(cross(c.p[m-1]-c.p[m-2], p[i]-c.p[m-2]))<=0) m--;
            c.p[m++] = p[i];
        }
        if (n > 1) m--;
    }
    int isPointIn(Point a) {
        p[n] = p[0];
        int w = 0;
        for (int i = 0; i < n; i++) {
            if (onSegment(a, p[i], p[i + 1])) {
                //if (p[i] == a || p[i + 1] == p) return 3;
                return 2;
            }
            int k = dcmp(cross(p[i + 1] - p[i], a - p[i]));
            int d1 = dcmp(p[i].y - a.y);
            int d2 = dcmp(p[i + 1].y - a.y);
            if(k > 0 && d1 <= 0 && d2 > 0) w++;
            if(k < 0 && d2 <= 0 && d1 > 0) w--;
        }
        return w != 0;
    }
    bool isConvex() {
        bool s[3] = {false, false, false};
        p[n] = p[0]; p[n + 1] = p[1];
        for (int i = 0; i < n; i++) {
             s[dcmp(cross(p[i + 1] - p[i], p[i + 2] - p[i])) + 1] = true;
             if (s[0] && s[2]) return false;
        }
        return true;
    }
    double cirum() {
        if (n == 1) return 0;
        //if (n == 2) return length(p[0] - p[1]);
        double sum = length(p[0] - p[n - 1]);
        for (int i = 1; i < n; i++)
            sum += length(p[i] - p[i - 1]);
        return sum;
    }
    double area() {
        double sum = 0;
        p[n] = p[0];
        for (int i = 0; i < n; i++) sum += cross(p[i], p[i + 1]);
        return sum * 0.5;
    }
    Point bary() {
        Point res(0, 0);
        double s = 0;
        p[n] = p[0];
        for (int i = 0; i < n; i++) {
            double t = cross(p[i], p[i + 1]);
            if (dcmp(t) == 0) continue;
            s += t;
            res = res + (p[i] + p[i + 1]) * (t / 3);
        }
        if (dcmp(s)) return res / s;
        return res;
    }
};

struct HLine {
    Point p, v;
    double ang;

    HLine() {}
    HLine(Point p, Point v) : p(p), v(v) {
        ang = atan2(v.y, v.x);
    }

    bool operator < (const HLine &u) const {
        return ang < u.ang;
    }
    bool onLeft(Point q) {
        return dcmp(cross(v, q - p)) > 0;
    }
    Point crossPoint(const HLine &u) {
        Point c = p - u.p;
        double t = cross(u.v, c) / cross(v, u.v);
        return p + v * t;
    }
};
struct Halfplanes {
    HLine hp[MP];
    int n;

    Point p[MP];
    int q[MP];
    int qh, qt;

    void clear() { n = 0; }
    void push(const HLine& h) {
        hp[n++] = h;
    }
    int halfplaneIntersect() {
        sort(hp, hp + n);
        qh = qt = 0;
        q[qt++] = 0;
        for (int i = 1; i < n; i++) {
            while (qt-qh > 1 && !hp[i].onLeft(p[qt-1])) qt--;
            while (qt-qh > 1 && !hp[i].onLeft(p[qh+1])) qh++;
            if (dcmp(cross(hp[i].v, hp[q[qt-1]].v)) == 0) {
                if(hp[q[qt-1]].onLeft(hp[i].p)) q[qt-1] = i;
            } else {
                q[qt++] = i;
            }
            if (qt-qh > 1) p[qt-1] = hp[q[qt-2]].crossPoint(hp[q[qt-1]]);
        }
        while (qh<qt && !hp[q[qh]].onLeft(p[qt - 1])) qt--;
        if (qt - qh <= 2) return 0;
        p[qh] = hp[q[qt - 1]].crossPoint(hp[q[qh]]);
        return qt - qh;
    }
    void getPolygon(Polygon& cv) {
        cv.n = 0;
        //if (!halfplaneIntersect()) return;
        for (int i = qh; i < qt; i++)
            cv.p[cv.n++] = p[i];
    }
};

struct Circle {
    Point c;
    double r;

    Circle() {}
    Circle(Point c, double r) : c(c), r(r) {}

    bool operator < (const Circle &C) const {
        return r < C.r;
    }
    bool operator == (const Circle &C) const {
        return c == C.c && dcmp(r - C.r) == 0;
    }

    Point pointAt(double a) {
        return Point(c.x + cos(a) * r, c.y + sin(a) * r);
    }
    int crossCircle(const Circle &C, Point &p1, Point &p2) {
        double d = length(c - C.c);
        if (dcmp(d) == 0) {
            if (dcmp(r - C.r) == 0) return -1; //两圆重合
            return 0;
        }
        if (dcmp(r + C.r - d) < 0) return 0;
        if (dcmp(fabs(r - C.r) - d) > 0) return 0;
        double a = atan2(C.c.y - c.y, C.c.x - c.x);
        double da = acos((r*r + d*d - C.r*C.r) / (2*r*d));
        p1 = pointAt(a - da);
        p2 = pointAt(a + da);
        if (dcmp(da) == 0) return 1;
        return 2;
    }
    int crossLine(Point a, Point b, Point &p1, Point &p2) {
        Point p = pointToLine(c, a, b);
        double d = length(p - c);
        if (dcmp(d - r) > 0) return 0;
        double t = sqrt(r*r - d*d) / length(a - b);
        Point v = (b - a) * t;
        p1 = p - v; p2 = p + v;
        if (dcmp(d - r) == 0) return 1;
        else return 2;
    }
    //切线方向向量
    int tangent(Point p, Point &v1, Point &v2) {
        Point u = c - p;
        double d = length(u);
        if (dcmp(d - r) < 0) return 0;
        if (dcmp(d - r) == 0) {
            v1 = v2 = rotate(u, pi * 0.5);
            return 1;
        }
        double ang = asin(r / d);
        v1 = rotate(u, -ang);
        v2 = rotate(u, ang);
        return 2;
        //p1 = p + v1 / v1.length() * sqrt(d*d - r*r); 切点
    }
    double areaCircle(Circle C) {
        double d = length(C.c - c);
        if (dcmp(r + C.r - d) <= 0) return 0;
        if (dcmp(r + d - C.r) <= 0) return pi * r * r;
        if (dcmp(C.r + d - r) <= 0) return pi * C.r * C.r;
        double d1 = (d*d + r*r - C.r*C.r) / (2*d);
        double d2 = (d*d + C.r*C.r - r*r) / (2*d);
        double h = sqrt(r*r - d1*d1);
        return acos(d1/r)*r*r - h*d1 + acos(d2/C.r)*C.r*C.r - h*d2;
    }
    //a到b的有向弧，a,b不一定在圆上
    double arc(Point a, Point b) {
        double cp = cross(a - c, b - c);
        double res=atan2(fabs(cp), dot(a - c, b - c));
        return cp > -eps ? res : -res;
    }
    //a到b的有向扇形面积
    double areaSector(Point a, Point b) {
        return arc(a, b) * r * r * 0.5;
    }
    //a到b的有向弓形面积
    double areaBow(Point a, Point b) {
        return areaSector(a, b) - cross(a - c, b - c) * 0.5;
    }
    //三角形另一个点是圆心，且圆心是原点
    double areaTriangle(Point a, Point b) {
        if(length(a) <= r && length(b) <= r) return cross(a, b) * 0.5;
        double d = disToLine(c, a, b);
        if(d >= r) return areaSector(a, b);
        Point p1, p2;
        crossLine(a, b, p1, p2);
        if(length(a) <= r) return cross(a, p2) * 0.5 + areaSector(p2, b);
        if(length(b) <= r) return areaSector(a, p1) + cross(p1, b) * 0.5;
        if(dcmp(dot(p1 - a, p1 - b)) >= 0) return areaSector(a, b);
        return areaSector(a, p1) + cross(p1, p2) * 0.5 + areaSector(p2, b);
    }
};
//外接圆
Circle cirsumCircle(Point A, Point B, Point C) {
    double bx = B.x - A.x, by = B.y - A.y;
    double cx = C.x - A.x, cy = C.y - A.y;
    double d = 2 * (bx*cy - by*cx);
    double x = (cy*(bx*bx+by*by) - by*(cx*cx+cy*cy)) / d + A.x;
    double y = (bx*(cx*cx+cy*cy) - cx*(bx*bx+by*by)) / d + A.y;
    Point p = Point(x, y);
    return Circle(p, length(A - p));
}
//内切圆
Circle inscribedCircle(Point A, Point B, Point C) {
    double a = length(B - C);
    double b = length(C - A);
    double c = length(A - B);
    Point p = (A*a + B*b + C*c) / (a+b+c);
    return Circle(p, disToLine(p, A, B));
}
