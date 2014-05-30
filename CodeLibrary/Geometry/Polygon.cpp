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
