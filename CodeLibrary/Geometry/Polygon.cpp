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

struct Circle {
    Point c;
    double r;

    Circle() {}
    Circle(Point c, double r) : c(c), r(r) {}

    bool operator < (const Circle &C) const {
        return r < C.r;
    }
    bool operator == (const Circle &C) const {
        return (c == C.c && dcmp(r - C.r) == 0);
    }

    Point pointAt(double a) {
        return Point(c.x + cos(a) * r, c.y + sin(a) * r);
    }
    int crossCircle(const Circle &C, Point &p1, Point &p2) {
        double d = (c - C.c).length();
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
    int crossLine(Line l, Point &p1, Point &p2) {
        Point p = l.pointToLine(c);
        double d = p.disToPoint(c);
        if (dcmp(d - r) > 0) return 0;
        double t = sqrt(r*r - d*d) / l.length();
        Point v = (l.b - l.a) * t;
        p1 = p - v; p2 = p + v;
        if (dcmp(d - r) == 0) return 1;
        else return 2;
    }
    //切线方向向量
    int tangent(Point p, Point &v1, Point &v2) {
        Point u = c - p;
        double d = u.length();
        if (dcmp(d - r) < 0) return 0;
        if (dcmp(d - r) == 0) {
            v1 = v2 = u.rotate(pi * 0.5);
            return 1;
        }
        double ang = asin(r / d);
        v1 = u.rotate(-ang);
        v2 = u.rotate(ang);
        return 2;
        //p1 = p + v1 / v1.length() * sqrt(d*d - r*r); 切点
    }
    double areaCircle(Circle C) {
        double d = (C.c - c).length();
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
        double cp=(a-c).cross(b-c);
        double res=atan2(fabs(cp), (a-c).dot(b-c));
        return cp > -eps ? res : -res;
    }
    //a到b的有向扇形面积
    double areaSector(Point a, Point b) {
        return arc(a, b) * r * r * 0.5;
    }
    //a到b的有向弓形面积
    double areaBow(Point a, Point b) {
        return areaSector(a, b) - (a-c).cross(b-c) * 0.5;
    }
    //三角形另一个点是圆心，且圆心是原点
    double areaTriangle(Point a, Point b) {
        if(a.length() <= r && b.length() <= r) return a.cross(b) * 0.5;
        Line ab(a, b);
        double d = ab.disToLine(c);
        if(d >= r) return areaSector(a, b);
        Point p1, p2;
        crossLine(ab, p1, p2);
        if(a.length() <= r) return a.cross(p2) * 0.5 + areaSector(p2, b);
        if(b.length() <= r) return areaSector(a,p1) + p1.cross(b) * 0.5;
        if(dcmp((p1-a).dot(p1-b)) >= 0) return areaSector(a, b);
        return areaSector(a, p1) + p1.cross(p2) * 0.5 + areaSector(p2, b);
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
    return Circle(p, (A - p).length());
}
//内切圆
Circle inscribedCircle(Point A, Point B, Point C) {
    double a = (B - C).length();
    double b = (C - A).length();
    double c = (A - B).length();
    Point p = (A*a + B*b + C*c) / (a+b+c);
    return Circle(p, Line(A, B).disToLine(p));
}

