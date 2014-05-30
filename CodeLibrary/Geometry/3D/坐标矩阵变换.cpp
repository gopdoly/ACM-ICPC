// ��P�����任Ϊ M[T-1] * ... * M[2] * M[1] * M[0] * P
// ���ݽ���ɣ��ȼ���mat = (M[T-1] * ... * M[0])�����P�任Ϊmat * P
// 4x4��α任����
struct Matrix {
    double v[4][4];

    // ����˷�
    inline Matrix operator * (const Matrix &b) const {
        Matrix c;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++) {
                c.v[i][j] = 0;
                for (int k = 0; k < 4; k++)
                    c.v[i][j] += v[i][k] * b.v[k][j];
            }
        return c;
    }

    inline Matrix pow(int k) {
        Matrix s, t;
        s.loadIdentity();
        t = *this;
        for (int i = k; i; i >>= 1) {
            if (i & 1) s = s * t;
            t = t * t;
        }
        return s;
    }

    // �任һ���㣬�൱���ҳ�������(x, y, z, 1}
    inline Point transform(Point P) const {
        double p[4] = {P.x, P.y, P.z, 1};
        double ans[4] = {0, 0, 0, 0};
        for(int i = 0; i < 4; i++)
            for(int k = 0; k < 4; k++)
                ans[i] += v[i][k] * p[k];
        return Point(ans[0], ans[1], ans[2]); // ans[3]�϶���1
    }

    // ��λ����
    void loadIdentity() {
        memset(v, 0, sizeof(v));
        v[0][0] = v[1][1] = v[2][2] = v[3][3] = 1;
    }

    // ƽ�ƾ���
    void loadTranslate(double a, double b, double c) {
        loadIdentity();
        v[0][3] = a; v[1][3] = b; v[2][3] = c;
    }

    // ���ž���
    void loadScale(double a, double b, double c) {
        loadIdentity();
        v[0][0] = a; v[1][1] = b; v[2][2] = c;
    }

    // �ƹ̶�����תһ���Ƕȵľ���
    // ��OP(P(a,b,c),P��O��)��ʱ��deg�Ƕ�
    void loadRotation(double a, double b, double c, double deg) {
        loadIdentity();
        double rad = deg / 180 * PI;
        double sine = sin(rad), cosine = cos(rad);
        double len = sqrt(a * a + b * b + c * c);
        Point L(a / len, b / len, c / len);
        v[0][0] = cosine + L.x * L.x * (1.0 - cosine);
        v[0][1] = L.x * L.y * (1 - cosine) - L.z * sine;
        v[0][2] = L.x * L.z * (1 - cosine) + L.y * sine;
        v[1][0] = L.y * L.x * (1 - cosine) + L.z * sine;
        v[1][1] = cosine + L.y * L.y * (1 - cosine);
        v[1][2] = L.y * L.z * (1 - cosine) - L.x * sine;
        v[2][0] = L.z * L.x * (1 - cosine) - L.y * sine;
        v[2][1] = L.z * L.y * (1 - cosine) + L.x * sine;
        v[2][2] = cosine + L.z * L.z * (1 - cosine);
    }
};

//uva12303 ��P413
const int maxn = 50000 + 10;
const int maxp = 50000 + 10;
Point3 P[maxn];
Plane planes[maxp];

int main() {
  int n, m, T;
  scanf("%d%d%d", &n, &m, &T);
  for(int i = 0; i < n; i++)
    scanf("%lf%lf%lf", &P[i].x, &P[i].y, &P[i].z);
  for(int i = 0; i < m; i++)
    scanf("%lf%lf%lf%lf", &planes[i].a, &planes[i].b, &planes[i].c, &planes[i].d);

  // ��P�����任Ϊ M[T-1] * ... * M[2] * M[1] * M[0] * P
  // ���ݽ���ɣ��ȼ���mat = (M[T-1] * ... * M[0])�����P�任Ϊmat * P
  Matrix4x4 mat;
  mat.loadIdentity();
  for(int i = 0; i < T; i++) {
    char op[100];
    double a, b, c, theta;
    scanf("%s%lf%lf%lf", op, &a, &b, &c);
    Matrix4x4 M;
    if(op[0] == 'T') M.loadTranslate(a, b, c);
    else if(op[0] == 'S') M.loadScale(a, b, c);
    else if(op[0] == 'R') { scanf("%lf", &theta); M.loadRotation(a, b, c, theta); }
    mat = M * mat;
  }

  // �任��
  for(int i = 0; i < n; i++) {
    Point3 ans = mat.transform(P[i]);
    printf("%.2lf %.2lf %.2lf\n", ans.x, ans.y, ans.z);
  }
  // �任ƽ��
  for(int i = 0; i < m; i++) {
    Point3 A[3];
    for(int j = 0; j < 3; j++) A[j] = mat.transform(planes[i].sample());
    Plane pl(A);
    printf("%.2lf %.2lf %.2lf %.2lf\n", pl.a, pl.b, pl.c, pl.d);
  }
  return 0;
}
