//向量A,B的夹角
double angle(Point a, Point b){
    return acos(dot(a, b) / length(a) / length(b));
}
//点p到平面p0-n(点法)的距离，n为单位向量时不用做除法
double disToPlane(Point p, Point p0, Point n) {
    return fabs(dot(p-p0, n)) / length(n);
}
//点p在平面p0-n上的投影，n为单位向量时不用做除法
Point getPlaneProjection(Point p, Point p0, Point n) {
    double d = dot(p-p0, n) / length(n);
    return p + n * d;
}
//直线p1-p2到平面p0-n的交点。假定交点唯一存在
Point linePlaneIntersection(Point p1, Point p2, Point p0, Point n) {
    Point v = p2 - p1;
    double t = dot(n, p0 - p1) / dot(n, p2 - p1); //判断分母是否为0，为0时平行
    return p1 + v * t; //如果是线段，判断t是不是在0和1之间
}
//三角形ABC面积的2倍
double area2(Point a, Point b, Point c) {
    return length(cross(b - a, c - a));
}
//点P在三角形ABC中
bool pointInTri(Point p, Point a, Point b, Point c) {
    double s = area2(a, b, c);
    double s1 = area2(p, a, b);
    double s2 = area2(p, a, c);
    double s3 = area2(p, b, c);
    return dcmp(s1 + s2 + s3 - s) == 0;
}
//三角形P0P1P2是否与线段AB相交
bool triSegIntersection(Point p0, Point p1, Point p2, Point a, Point b, Point& p) {
    Point n = cross(p1 - p0, p2 - p0);
    if (dcmp(dot(n, b - a)) == 0) return false;
    double t = dot(n, p0 - a) / dot(n, b - a);
    if (dcmp(t) < 0 || dcmp(t - 1) > 0) return false;
    p = a + (b - a) * t;
    return pointInTri(p, p0, p1, p2);
}
//点P到直线AB的距离
double disToLine(Point p, Point a, Point b) {
    Point v1 = b - a, v2 = p - a;
    return length(cross(v1, v2)) / length(v1);
}
//点P到线段AB的距离
double disToSeg(Point p, Point a, Point b) {
    if (a == b) return length(p - a);
    Point v1 = b - a, v2 = p - a, v3 = p - b;
    if (dcmp(dot(v1, v2)) < 0) return length(v2);
    else if (dcmp(dot(v1, v3)) > 0) return length(v3);
    else return length(cross(v1, v2)) / length(v1);
}
//返回AB, AC, AD呈右手系的混合积。也等于四面体ABCD的有向体积的6倍
double volume6(Point a, Point b, Point c, Point d) {
    return dot(d-a, cross(b-a, c-a));
}
//求异面直线p1+su和p2+tv的公垂线对应的s。如果平行/重合，返回false
bool lineDistance(Point p1, Point u, Point p2, Point v, double& s) {
    double b = dot(u, u) * dot(v, v) - dot(u, v) * dot(u, v);
    if (dcmp(b) == 0) return false;
    double a = dot(u, v) * dot(v, p1 - p2) - dot(v, v) * dot(u, p1 - p2);
    s = a / b;
    return true;
}
//两条直线的距离
double lineDisToLine(Point p1, Point v1, Point p2, Point v2) {
    Point u = cross(v1, v2), a = p1 - p2;
    if (dcmp(length(u)) == 0) return length(cross(a, v2)) / length(v2);
    return fabs(dot(a, u)) / length(u);
}
//点A绕OP轴逆时针(O向P看)旋转ang度
Point rotate(Point a, Point p, double ang) {
    p = p / length(p);
    Point v = cross(a, p);
    double d = dot(a, p);
    double sn = sin(ang), cs = cos(ang);
    return a * cs - v * sn + p * d * (1 - cs);
}
//平面O-n1上的点a，旋转到平面O-n2上(O为原点(0,0,0))
Point rotate2(Point a, Point n1, Point n2) {
    Point u = cross(n1, n2);
    if (dcmp(length(u)) == 0) return a;
    double rad = acos(dot(n1, n2) / length(n1) / length(n2));
    return rotate(a, u, rad);
}
