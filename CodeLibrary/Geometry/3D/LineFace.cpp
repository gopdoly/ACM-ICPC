//����A,B�ļн�
double angle(Point a, Point b){
    return acos(dot(a, b) / length(a) / length(b));
}
//��p��ƽ��p0-n(�㷨)�ľ��룬nΪ��λ����ʱ����������
double disToPlane(Point p, Point p0, Point n) {
    return fabs(dot(p-p0, n)) / length(n);
}
//��p��ƽ��p0-n�ϵ�ͶӰ��nΪ��λ����ʱ����������
Point getPlaneProjection(Point p, Point p0, Point n) {
    double d = dot(p-p0, n) / length(n);
    return p + n * d;
}
//ֱ��p1-p2��ƽ��p0-n�Ľ��㡣�ٶ�����Ψһ����
Point linePlaneIntersection(Point p1, Point p2, Point p0, Point n) {
    Point v = p2 - p1;
    double t = dot(n, p0 - p1) / dot(n, p2 - p1); //�жϷ�ĸ�Ƿ�Ϊ0��Ϊ0ʱƽ��
    return p1 + v * t; //������߶Σ��ж�t�ǲ�����0��1֮��
}
//������ABC�����2��
double area2(Point a, Point b, Point c) {
    return length(cross(b - a, c - a));
}
//��P��������ABC��
bool pointInTri(Point p, Point a, Point b, Point c) {
    double s = area2(a, b, c);
    double s1 = area2(p, a, b);
    double s2 = area2(p, a, c);
    double s3 = area2(p, b, c);
    return dcmp(s1 + s2 + s3 - s) == 0;
}
//������P0P1P2�Ƿ����߶�AB�ཻ
bool triSegIntersection(Point p0, Point p1, Point p2, Point a, Point b, Point& p) {
    Point n = cross(p1 - p0, p2 - p0);
    if (dcmp(dot(n, b - a)) == 0) return false;
    double t = dot(n, p0 - a) / dot(n, b - a);
    if (dcmp(t) < 0 || dcmp(t - 1) > 0) return false;
    p = a + (b - a) * t;
    return pointInTri(p, p0, p1, p2);
}
//��P��ֱ��AB�ľ���
double disToLine(Point p, Point a, Point b) {
    Point v1 = b - a, v2 = p - a;
    return length(cross(v1, v2)) / length(v1);
}
//��P���߶�AB�ľ���
double disToSeg(Point p, Point a, Point b) {
    if (a == b) return length(p - a);
    Point v1 = b - a, v2 = p - a, v3 = p - b;
    if (dcmp(dot(v1, v2)) < 0) return length(v2);
    else if (dcmp(dot(v1, v3)) > 0) return length(v3);
    else return length(cross(v1, v2)) / length(v1);
}
//����AB, AC, AD������ϵ�Ļ�ϻ���Ҳ����������ABCD�����������6��
double volume6(Point a, Point b, Point c, Point d) {
    return dot(d-a, cross(b-a, c-a));
}
//������ֱ��p1+su��p2+tv�Ĺ����߶�Ӧ��s�����ƽ��/�غϣ�����false
bool lineDistance(Point p1, Point u, Point p2, Point v, double& s) {
    double b = dot(u, u) * dot(v, v) - dot(u, v) * dot(u, v);
    if (dcmp(b) == 0) return false;
    double a = dot(u, v) * dot(v, p1 - p2) - dot(v, v) * dot(u, p1 - p2);
    s = a / b;
    return true;
}
//����ֱ�ߵľ���
double lineDisToLine(Point p1, Point v1, Point p2, Point v2) {
    Point u = cross(v1, v2), a = p1 - p2;
    if (dcmp(length(u)) == 0) return length(cross(a, v2)) / length(v2);
    return fabs(dot(a, u)) / length(u);
}
//��A��OP����ʱ��(O��P��)��תang��
Point rotate(Point a, Point p, double ang) {
    p = p / length(p);
    Point v = cross(a, p);
    double d = dot(a, p);
    double sn = sin(ang), cs = cos(ang);
    return a * cs - v * sn + p * d * (1 - cs);
}
//ƽ��O-n1�ϵĵ�a����ת��ƽ��O-n2��(OΪԭ��(0,0,0))
Point rotate2(Point a, Point n1, Point n2) {
    Point u = cross(n1, n2);
    if (dcmp(length(u)) == 0) return a;
    double rad = acos(dot(n1, n2) / length(n1) / length(n2));
    return rotate(a, u, rad);
}
