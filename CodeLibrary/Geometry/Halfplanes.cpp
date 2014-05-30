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
