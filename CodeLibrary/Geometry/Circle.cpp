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
