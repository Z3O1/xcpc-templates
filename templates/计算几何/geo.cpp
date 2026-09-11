
namespace Geo {

const db eps = 1e-10;

struct p2;
struct seg;
int sign(db x) { return x < -eps ? -1 : x > eps; }
int cmp(db x, db y) { return sign(x - y); }
const db pi = acos(db(-1));
struct p2 {
    db x, y;
    db det(const p2 &b) const { return x * b.y - y * b.x; }
    db alpha() const { return atan2(y, x); }
    bool operator==(const p2 &b) const {
        return !cmp(x, b.x) && !cmp(y, b.y);
    }
    bool operator<(const p2 &b) const {
        int c = cmp(x, b.x);
        if(c) return !~c;
        return !~cmp(y, b.y);
    }
};
p2 r90(p2 x) { return {-x.y, x.x}; }
p2 rot(p2 x, db a) { return {x.x * cos(a) - x.y * sin(a), x.x * sin(a) + x.y * cos(a)}; }
p2 operator+(p2 x, p2 y) { return {x.x + y.x, x.y + y.y}; }
p2 operator-(p2 x, p2 y) { return {x.x - y.x, x.y - y.y}; }
p2 operator/(p2 x, db y) { return {x.x / y, x.y / y}; }
p2 operator*(p2 x, db y) { return {x.x * y, x.y * y}; }
p2 operator*(db y, p2 x) { return {x.x * y, x.y * y}; }
db operator*(p2 x, p2 y) { return x.x * y.x + x.y * y.y; }
p2 operator+=(p2 &x, p2 y) { return x = x + y; }
p2 operator-=(p2 &x, p2 y) { return x = x + y; }
p2 operator/=(p2 &x, p2 y) { return x = x + y; }
db dis2(p2 x) { return x.x * x.x + x.y * x.y; }
db dis(p2 x) { return hypot(x.x, x.y); }
p2 unit(p2 x) { return x / dis(x); }
struct seg {
    p2 x, y;
    p2 dir() const { return y - x; }
};
db cross(const seg &p, const p2 &q) { return p.dir().det(q - p.x); }
int crossop(const seg &p, const p2 &q) { return sign(cross(p, q)); }
bool inc(seg l, const p2 &b) { return !crossop(l, b); }
bool eqll(const seg &a, const seg &b) {
    return inc(a, b.x) && inc(a, b.y);
}
bool chkll(const seg &a, const seg &b) {
    db a1 = cross(b, a.x), a2 = -cross(b, a.y);
    return !!sign(a1 + a2);
}
p2 isll(const seg &a, const seg &b) {
    db a1 = cross(b, a.x), a2 = -cross(b, a.y);
    return (a.x * a2 + a.y * a1) / (a1 + a2);
}
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
bool chkss_s(const seg &p, const seg &q) {
    p2 p1 = p.x, p2 = p.y, q1 = q.x, q2 = q.y;
    return crossop(p, q1) * crossop(p, q2) < 0 && crossop(q, p1) * crossop(q, p2) < 0;
}
bool isMid(const db &a, const db &m, const db &b) { return !cmp(a, b) || !cmp(b, m) || cmp(a, m) != cmp(b, m); }
bool isMid(const p2 &a, const p2 &m, const p2 &b) { return isMid(a.x, m.x, b.x) && isMid(a.y, m.y, b.y); }
bool ons(const seg &p, const p2 &q) { return !crossop(p, q) && isMid(p.x, q, p.y); }
bool ons_s(const seg &p, const p2 &q) {
    return !crossop(p, q) && sign((q - p.x) * (p.x - p.y)) * sign((q - p.y) * (p.x - p.y)) < 0;
}
p2 proj(const seg &p, const p2 &q) {
    p2 dir = p.y - p.x;
    return p.x + dir * ((dir * (q - p.x)) / dis2(dir));
}
p2 reflect(const seg &p, const p2 &q) {
    return proj(p, q) * 2 - q;
}
db nearest(const seg &p, const p2 &q) {
    p2 h = proj(p, q);
    if(isMid(p.x, h, p.y)) return dis2(q - h);
    return min(dis2(p.x - q), dis2(p.y - q));
}
db disss(const seg &p, const seg &q) {
    if(chkss(p, q)) return 0;
    return min({nearest(p, q.x), nearest(p, q.y), nearest(q, p.x), nearest(q, p.y)});
}
db area(int n, p2 *a) {
    db ret = 0;
    For(i, 0, n - 1) ret += a[i].det(a[(i + 1) % n]);
    return ret / 2;
}
int contain(int n, p2 *a, const p2 &p) {  // 0: outside, 1: on seg, 2: inside
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
db convex_diameter(int n, p2 *a) {
    if(n <= 1) return 0;
    int is = 0, js = 0;
    For(k, 1, n - 1) is = a[k] < a[is] ? k : is, js = a[js] < a[k] ? k : js;
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
vector<p2> convex_cut(int n, p2 *a, const seg &q) {
    vector<p2> qs;
    For(i, 0, n - 1) {
        p2 p1 = a[i], p2 = a[(i + 1) % n];
        int d1 = crossop(q, p1), d2 = crossop(q, p2);
        if(d1 >= 0) qs.push_back(p1);
        if(d1 * d2 < 0) qs.push_back(isll({p1, p2}, q));
    }
    return qs;
}
db circumcircle_diameter(p2 x, p2 y, p2 z) {
    if(!crossop({x, y}, z)) return -1;
    return dis(x - y) * dis(y - z) * dis(z - x) / abs((y - x).det(z - x));
}
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
