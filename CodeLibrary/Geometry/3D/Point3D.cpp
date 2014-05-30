struct Point {
    double x, y, z;

    Point() {}
    Point(double x, double y, double z)
    : x(x), y(y), z(z) { }

    Point operator + (const Point &a) const {
        return Point(x + a.x, y + a.y, z + a.z);
    }
    Point operator - (const Point &a) const {
        return Point(x - a.x, y - a.y, z - a.z);
    }
    Point operator * (double k) const {
        return Point(x * k, y * k, z * k);
    }
    Point operator / (double k) const {
        return Point(x / k, y / k, z / k);
    }
    void read() {
        scanf("%lf%lf%lf", &x, &y, &z);
    }
    void show() {
        printf("%f %f %f\n", x, y, z);
    }
};
double dot(const Point& a, const Point& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
Point cross(const Point& a, const Point& b) {
    return Point(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
}
double length(const Point& a) {
    return sqrt(dot(a, a));
}
