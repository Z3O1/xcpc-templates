// Geo(): 计算几何 · 基础几何库 (点/线段/多边形/圆)

// eps = 1e-10 是**绝对**容差:sign/cmp 与所有位置/相交判定都走它;坐标量级很大(>1e6)时请自行放宽。
// p2 兼作点与向量;seg{x, y} 是有向线段(方向 x -> y),也可当直线用。
namespace Geo {

const db eps = 1e-10;

struct p2;
struct seg;
// 符号(带 eps 容差):x < -eps 为 -1,x > eps 为 1,否则 0。
int sign(db x) { return x < -eps ? -1 : x > eps; }
// x 与 y 的容差比较:返回 x-y 的符号。
int cmp(db x, db y) { return sign(x - y); }
const db pi = acos(db(-1));
struct p2 {
    db x, y;
    // 叉积 x×b(带号面积 ×2)。
    db det(const p2 &b) const { return x * b.y - y * b.x; }
    // 极角,取值 (-pi, pi]。
    db alpha() const { return atan2(y, x); }
    // 按 eps 容差判等。
    bool operator==(const p2 &b) const {
        return !cmp(x, b.x) && !cmp(y, b.y);
    }
    // 字典序(x 优先),用于排序/二分。
    bool operator<(const p2 &b) const {
        int c = cmp(x, b.x);
        if(c) return c < 0;
        return cmp(y, b.y) < 0;
    }
};
// 逆时针旋转 90 度。
p2 r90(p2 x) { return {-x.y, x.x}; }
// 绕原点逆时针旋转 a 弧度。
p2 rot(p2 x, db a) { return {x.x * cos(a) - x.y * sin(a), x.x * sin(a) + x.y * cos(a)}; }
// 向量相加。
p2 operator+(p2 x, p2 y) { return {x.x + y.x, x.y + y.y}; }
// 向量相减。
p2 operator-(p2 x, p2 y) { return {x.x - y.x, x.y - y.y}; }
// 向量除以标量;要求 y 非零。
p2 operator/(p2 x, db y) { return {x.x / y, x.y / y}; }
// 数乘(标量在右)。
p2 operator*(p2 x, db y) { return {x.x * y, x.y * y}; }
// 数乘(标量在左)。
p2 operator*(db y, p2 x) { return {x.x * y, x.y * y}; }
// 点积。
db operator*(p2 x, p2 y) { return x.x * y.x + x.y * y.y; }
// 向量自加,返回新值。
p2 operator+=(p2 &x, p2 y) { return x = x + y; }
// 向量自减,返回新值。
p2 operator-=(p2 &x, p2 y) { return x = x - y; }
// 向量自除,返回新值;要求 y 非零。
p2 operator/=(p2 &x, db y) { return x = x / y; }
// 模长平方。
db dis2(p2 x) { return x.x * x.x + x.y * x.y; }
// 模长。
db dis(p2 x) { return hypot(x.x, x.y); }
// 单位化;要求 x 非零(零向量得 inf/nan)。
p2 unit(p2 x) { return x / dis(x); }
struct seg {
    p2 x, y;
    // 方向向量 y - x。
    p2 dir() const { return y - x; }
};
// 有向直线 p 与点 q 的叉积(q 在左 > 0)。
db cross(const seg &p, const p2 &q) { return p.dir().det(q - p.x); }
// 上述叉积的符号。
int crossop(const seg &p, const p2 &q) { return sign(cross(p, q)); }
// 点 b 是否在有向直线 l 上。
bool inc(seg l, const p2 &b) { return !crossop(l, b); }
// 两直线是否重合(方向可反)。
bool eqll(const seg &a, const seg &b) {
    return inc(a, b.x) && inc(a, b.y);
}
// 两直线是否不平行(有唯一交点)。
bool chkll(const seg &a, const seg &b) {
    db a1 = cross(b, a.x), a2 = -cross(b, a.y);
    return !!sign(a1 + a2);
}
// 两直线交点;要求不平行(见 chkll)。
p2 isll(const seg &a, const seg &b) {
    db a1 = cross(b, a.x), a2 = -cross(b, a.y);
    return (a.x * a2 + a.y * a1) / (a1 + a2);
}
// 两线段是否相交(含端点相接、共线重叠)。
bool chkss(const seg &p, const seg &q) {
    auto is = [](db l1, db r1, db l2, db r2) {
        if(l1 > r1) swap(l1, r1);
        if(l2 > r2) swap(l2, r2);
        return !(cmp(r1, l2) == -1 || cmp(r2, l1) == -1);
    };
    p2 p1 = p.x, p2 = p.y, q1 = q.x, q2 = q.y;
    return is(p1.x, p2.x, q1.x, q2.x) && is(p1.y, p2.y, q1.y, q2.y) &&
           crossop(p, q1) * crossop(p, q2) <= 0 && crossop(q, p1) * crossop(q, p2) <= 0;
}
// 两线段是否严格相交(交点严格在两段内部,端点/共线不算)。
bool chkss_s(const seg &p, const seg &q) {
    p2 p1 = p.x, p2 = p.y, q1 = q.x, q2 = q.y;
    return crossop(p, q1) * crossop(p, q2) < 0 && crossop(q, p1) * crossop(q, p2) < 0;
}
// m 是否落在 a、b 闭区间内(eps 容差)。
bool isMid(const db &a, const db &m, const db &b) { return cmp(a, m) * cmp(b, m) <= 0; }
// 逐分量判定 m 是否落在 a、b 之间。
bool isMid(const p2 &a, const p2 &m, const p2 &b) { return isMid(a.x, m.x, b.x) && isMid(a.y, m.y, b.y); }
// 点 q 是否在线段 p 上(含端点)。
bool ons(const seg &p, const p2 &q) { return !crossop(p, q) && isMid(p.x, q, p.y); }
// 点 q 是否严格在线段 p 内部(不含端点)。
bool ons_s(const seg &p, const p2 &q) {
    return !crossop(p, q) && sign((q - p.x) * (p.x - p.y)) * sign((q - p.y) * (p.x - p.y)) < 0;
}
// 点 q 在直线 p 上的投影;要求直线非零长。
p2 proj(const seg &p, const p2 &q) {
    p2 dir = p.y - p.x;
    return p.x + dir * ((dir * (q - p.x)) / dis2(dir));
}
// 点 q 关于直线 p 的对称点;要求直线非零长。
p2 reflect(const seg &p, const p2 &q) {
    return proj(p, q) * 2 - q;
}
// 点 q 到线段 p 的最近距离平方;线段可退化为点。
db nearest(const seg &p, const p2 &q) {
    if(dis2(p.dir()) == 0) return dis2(q - p.x);
    p2 h = proj(p, q);
    if(isMid(p.x, h, p.y)) return dis2(q - h);
    return min(dis2(p.x - q), dis2(p.y - q));
}
// 两线段间最近距离(平方);相交(含端点/共线)返回 0,支持退化线段。
db disss(const seg &p, const seg &q) {
    if(chkss(p, q)) return 0;
    return min({nearest(p, q.x), nearest(p, q.y), nearest(q, p.x), nearest(q, p.y)});
}
// 多边形带号面积(顶点顺序任意,顺/逆时针符号相反)。
db area(int n, p2 *a) {
    db ret = 0;
    For(i, 0, n - 1) ret += a[i].det(a[(i + 1) % n]);
    return ret / 2;
}
// 点 p 与多边形 a 的关系(射线法/奇偶规则):0 严格在外,1 在边界,2 严格在内;支持凹多边形。
int contain(int n, p2 *a, const p2 &p) {
    int ret = 0;
    For(i, 0, n - 1) {
        p2 u = a[i], v = a[(i + 1) % n];
        if(ons({u, v}, p)) return 1;
        if(cmp(u.y, v.y) <= 0) swap(u, v);
        if(cmp(p.y, u.y) > 0 || cmp(p.y, v.y) <= 0) continue;
        ret ^= crossop({p, u}, v) > 0;
    }
    return ret * 2;
}
// 求凸包,返回点数,k 个顶点按逆时针存在 b[0..k-1]。nos=0 去掉共线点,nos=1 保留。
// **前置条件:输出缓冲 b 至少要能放 2n + 2 个点**(不是 n 个!扫下凸壳时会先写满 n 个位置,
// 之后上凸壳还会再写)。
// 输入 a 会被原地排序;b 与 a 不能是同一块内存。
int convex_hull(int n, p2 *a, p2 *b, bool nos = 0) {
    if(n <= 1) return 1;
    sort(a, a + n);
    int k = 0;
    For(i, 0, n - 1) {
        while(k > 1 && crossop({b[k - 2], b[k - 1]}, a[i]) <= -nos) --k;
        b[k++] = a[i];
    }
    int t = k;
    rFor(i, n - 2, 0) {
        while(k > t && crossop({b[k - 2], b[k - 1]}, a[i]) <= -nos) --k;
        if(i) b[k++] = a[i];
    }
    return k;
}
// 凸多边形直径(最远点对)的距离平方;要求 a 逆时针、n >= 3(旋转卡壳)。
db convex_diameter(int n, p2 *a) {
    if(n <= 1) return 0;
    int is = 0, js = 0;
    For(k, 1, n - 1) {
        if(a[k] < a[is]) is = k;
        if(a[js] < a[k]) js = k;
    }
    int i = is, j = js;
    db ret = dis2(a[i] - a[j]);
    do {
        if((a[(i + 1) % n] - a[i]).det(a[(j + 1) % n] - a[j]) >= 0)
            (++j) %= n;
        else
            (++i) %= n;
        ret = max(ret, dis2(a[i] - a[j]));
    } while(i != is || j != js);
    return ret;
}
// 用直线 q 裁凸多边形,保留 q 左侧(cross(q, ·) >= 0);返回裁后多边形。
vector<p2> convex_cut(int n, p2 *a, const seg &q) {
    vector<p2> qs;
    For(i, 0, n - 1) {
        p2 u = a[i], v = a[(i + 1) % n];
        int d1 = crossop(q, u), d2 = crossop(q, v);
        if(d1 >= 0) qs.push_back(u);
        if(d1 * d2 < 0) qs.push_back(isll({u, v}, q));
    }
    return qs;
}
// 三点外接圆直径;三点共线返回 -1。
db circumcircle_diameter(p2 x, p2 y, p2 z) {
    if(!crossop({x, y}, z)) return -1;
    return dis(x - y) * dis(y - z) * dis(z - x) / abs((y - x).det(z - x));
}
// 三点外接圆圆心;要求三点不共线(否则除零)。
array<db, 2> circumcenter(p2 x, p2 y, p2 z) {
    auto [x1, y1] = x;
    auto [x2, y2] = y;
    auto [x3, y3] = z;
    auto A = dis2(x), B = dis2(y), C = dis2(z);
    auto u1 = x1 - x2, u2 = x1 - x3, v1 = y1 - y2, v2 = y1 - y3;
    return {db((C - A) * v1 - (B - A) * v2) / (2 * u1 * v2 - 2 * u2 * v1), db((C - A) * u1 - (B - A) * u2) / (2 * v1 * u2 - 2 * v2 * u1)};
}

}  // namespace Geo
using namespace Geo;
