//三维凸包p290
struct Face{
    int v[3];
    Face(int a, int b, int c) {
        v[0] = a; v[1] = b; v[2] = c;
    }
    Point normal(Point* p) const {
        return cross(p[v[1]]-p[v[0]], p[v[2]]-p[v[0]]);
    }
    bool cansee(Point* p, int i) const {
        return dot(p[i]-p[v[0]], normal(p)) > 0;
    }
};
//扰动 #include <cstdlib>
double rand01() { return rand() / (double)RAND_MAX; } //0-1的随机数
double randeps() {return (rand01() - 0.5) * eps; } //-eps/2到eps/2的随机数
Point add_noise(const Point& p) {
    return Point(p.x+randeps(), p.y+randeps(), p.z+randeps());
}
//增量法求三维凸包。
//没有考虑各种特殊情况(如4点共面)。实践中，调用前对输入点进行微小扰动
bool vis[N][N];
vector<Face> CH3D (Point* p, int n) {
    memset(vis, 0, sizeof(vis));
    vector<Face> cur;
    cur.push_back(Face(0, 1, 2));
    cur.push_back(Face(2, 1, 0));
    for (int i = 3; i < n; i++) {
        vector<Face> next;
        int sz = cur.size();
        for (int j = 0; j < sz; j++) {
            Face& f = cur[j];
            bool res = f.cansee(p, i);
            if (!res) next.push_back(f);
            for (int k = 0; k < 3; k++) vis[f.v[k]][f.v[(k+1)%3]] = res;
        }
        for (int j = 0; j < sz; j++) {
            for (int k = 0; k < 3; k++) {
                int a = cur[j].v[k], b = cur[j].v[(k+1)%3];
                if (vis[a][b] != vis[b][a] && vis[a][b])
                    next.push_back(Face(a, b, i));
            }
        }
        cur=next;
    }
    return cur;
}

//getFace()求三维凸包
//getBary()求重心
struct Polyhedron {
    Point p[N], tp[N];
    int n;
    vector<Face> fc;

    void read(int np) {
        n = np;
        for (int i = 0; i < n; i++) p[i] = read_point();
    }
    void getFace() {
        //输入点不判重，扰动后很容易出现极小的面
        for (int i = 0; i < n; i++) tp[i] = add_noise(p[i]);
        fc = CH3D(tp, n);
    }
    Point getBary() {
        Point res;
        double vol = 0;
        int sz = fc.size();
        for (int i = 0; i < sz; i++) {
            int a = fc[i].v[0], b = fc[i].v[1], c = fc[i].v[2];
            double tmp = dot(p[c], cross(p[a], p[b]));
            vol += tmp;
            res = res + (p[a] + p[b] + p[c]) * (tmp * 0.25);
        }
        res = res / vol;
        return res;
    }
};