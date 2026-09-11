// geo 自测:p2 点/向量运算、比较函数、距离、线段位置关系、凸包、外心等与手算/独立参考对照。
//
// ⚠ 本 check 目前是 FAIL 状态(预期),原因是 templates/计算几何/geo.cpp 第 33、34 行:
//     p2 operator-=(p2 &x, p2 y) { return x = x + y; }   // 应为 x = x - y
//     p2 operator/=(p2 &x, p2 y) { return x = x + y; }   // 应为 x = x / y,且第二参数应为 db
//   也就是「减等」与「除等」都写成了加法。本文件末尾「已知 bug 复现」段会先把这两个运算符的
//   实测值原样打印出来,再断言失败(exit 1)。修好 geo.cpp 第 33/34 行后本文件应变成 PASSED;
//   若仓库主人选择直接删掉 /= 重载,需要同步改本文件最后一段。
//
//   注:第 34 行的签名只有 (p2&, p2),连 (p2&, db) 版都没有 —— 写 `a /= 2.0` 会直接编译失败
//   (error: no match for 'operator/='),所以这里用 SFINAE 探测该重载是否存在,再单独测
//   (p2&, p2) 版的实测行为,以保证 check 自身在修 bug 前后都能编译通过。
//
// 本次运行结论:+= 正确;其他接口(见下面各段 [ok])全部与手算或独立参考实现一致;
//   只有 -= 与 /= 复现了上面的已知 bug,故本文件在 bug 修好前恒为 FAIL。
#include "../_check_base.hpp"
#include "geo.cpp"

// ---- 小工具 ----
static p2 P(db x, db y) { return {x, y}; }
static db Abs(db x) { return x < 0 ? -x : x; }
static bool eqd(db a, db b, db tol = 1e-9) { return Abs(a - b) < tol; }
static bool eqp(p2 a, p2 b, db tol = 1e-9) { return eqd(a.x, b.x, tol) && eqd(a.y, b.y, tol); }
static db dd(p2 a, p2 b) { return dis2(a - b); }
static void pr(p2 a) { printf("(%g,%g)", (double) a.x, (double) a.y); }

// 探测 /= 的两个候选重载是否存在(不存在的那个不能直接写,否则整个 check 编译不过)
// 注意:表达式里必须写 T 而不是直接写 p2 —— 否则它不依赖模板参数、不走 SFINAE,会直接编译报错
template <class T, class = void> struct has_diveq_db : std::false_type {};
template <class T> struct has_diveq_db<T, std::void_t<decltype(std::declval<T &>() /= std::declval<db>())>> : std::true_type {};
template <class T, class = void> struct has_diveq_p2 : std::false_type {};
template <class T> struct has_diveq_p2<T, std::void_t<decltype(std::declval<T &>() /= std::declval<p2>())>> : std::true_type {};
// 真正调用也必须放在模板里,而且操作数类型要依赖模板参数:否则表达式在模板定义时就被解析,
// 被 if constexpr 丢弃的分支照样报 "no match for operator/="
template <class T = p2> bool run_diveq_db(T a, db s, T &out) {
    if constexpr (has_diveq_db<T>::value) {
        out = a, out /= s;
        return true;
    } else {
        (void) a, (void) s, (void) out;
        return false;
    }
}
template <class T = p2, class Y = p2> bool run_diveq_p2(T a, Y y, T &out) {  // 返回该重载是否存在
    if constexpr (has_diveq_p2<T>::value) {
        out = a, out /= y;
        return true;
    } else {
        (void) a, (void) y, (void) out;
        return false;
    }
}

// ---- 独立参考实现(整数坐标下精确) ----
static db refNearest(p2 a, p2 b, p2 q) {  // 点到线段距离平方
    p2 dir = b - a;
    db t = (dir * (q - a)) / dis2(dir);
    return dd(q, a + dir * max((db) 0, min((db) 1, t)));
}
static db refDisss(p2 a, p2 b, p2 c, p2 e) {  // 线段间距离平方(相交记 0)
    if(chkss(seg{a, b}, seg{c, e})) return 0;
    return min({refNearest(a, b, c), refNearest(a, b, e), refNearest(c, e, a), refNearest(c, e, b)});
}
static bool refChkss(p2 a, p2 b, p2 c, p2 e) {  // 闭线段相交(参数法,行列式精确)
    p2 d = b - a, f = e - c, g = c - a;
    db den = d.det(f);
    if(Abs(den) > 1e-12) {
        db t = g.det(f) / den, u = g.det(d) / den;
        return t >= -1e-12 && t <= 1 + 1e-12 && u >= -1e-12 && u <= 1 + 1e-12;
    }
    if(Abs(g.det(d)) > 1e-12) return false;  // 平行且不共线
    db l1 = min(a.x, b.x), r1 = max(a.x, b.x), l2 = min(c.x, e.x), r2 = max(c.x, e.x);
    if(Abs(a.x - b.x) < 1e-12) l1 = min(a.y, b.y), r1 = max(a.y, b.y), l2 = min(c.y, e.y), r2 = max(c.y, e.y);
    return !(r1 < l2 - 1e-12 || r2 < l1 - 1e-12);
}
static vector<p2> refHull(vector<p2> v) {  // Andrew 单调链(去重、去共线)
    sort(v.begin(), v.end(), [](p2 a, p2 b) { return a.x != b.x ? a.x < b.x : a.y < b.y; });
    v.erase(unique(v.begin(), v.end(), [](p2 a, p2 b) { return a.x == b.x && a.y == b.y; }), v.end());
    int n = v.size();
    if(n <= 1) return v;
    vector<p2> h(2 * n);
    int k = 0;
    ForD(i, 0, n) {
        while(k > 1 && (h[k - 1] - h[k - 2]).det(v[i] - h[k - 2]) <= 0) --k;
        h[k++] = v[i];
    }
    int t = k;
    rFor(i, n - 2, 0) {
        while(k > t && (h[k - 1] - h[k - 2]).det(v[i] - h[k - 2]) <= 0) --k;
        h[k++] = v[i];
    }
    h.resize(k - 1);
    return h;
}
static int refContain(const vector<p2> &h, p2 q) {  // 逆时针凸多边形:0 外 / 1 边界 / 2 内
    int out = 0, on = 0;
    ForD(i, 0, h.size()) {
        p2 u = h[i], v = h[(i + 1) % h.size()];
        db c = (v - u).det(q - u);
        if(c == 0 && (q - u) * (q - v) <= 0) on = 1;
        if(c < 0) out = 1;
    }
    return on ? 1 : (out ? 0 : 2);
}
static vector<p2> refCut(const vector<p2> &h, seg q) {  // 半平面裁剪(保留 cross >= 0 一侧)
    vector<p2> r;
    ForD(i, 0, h.size()) {
        p2 p1 = h[i], p2r = h[(i + 1) % h.size()];
        db c1 = cross(q, p1), c2 = cross(q, p2r);
        if(c1 >= -1e-12) r.push_back(p1);
        if((c1 > 1e-12 && c2 < -1e-12) || (c1 < -1e-12 && c2 > 1e-12)) r.push_back(p1 + (p2r - p1) * (c1 / (c1 - c2)));
    }
    return r;
}
static db polyArea(const vector<p2> &h) {  // 无符号面积
    db s = 0;
    ForD(i, 0, h.size()) s += h[i].det(h[(i + 1) % h.size()]);
    return Abs(s) / 2;
}

int main() {
    // ===== 1. eps / pi / sign / cmp =====
    CHECK(eps == 1e-10, "eps == 1e-10");
    CHECK(eqd(pi, acosl(-1.0L)) && eqd(pi, 3.14159265358979L), "pi == acos(-1)");
    CHECK(sign(0) == 0 && sign(eps / 2) == 0 && sign(-eps / 2) == 0 && sign(eps) == 0, "sign:|x| <= eps 归零");
    CHECK(sign(eps * 10) == 1 && sign(-eps * 10) == -1, "sign:刚超出 eps 取符号");
    CHECK(sign(1e-6L) == 1 && sign(-1e-6L) == -1 && sign(100) == 1 && sign(-100) == -1, "sign:大值取符号");
    CHECK(cmp(1, 1) == 0 && cmp(1 + eps / 4, 1) == 0 && cmp(1, 1 + eps / 4) == 0, "cmp:容差内为 0");
    CHECK(cmp(2, 1) == 1 && cmp(1, 2) == -1 && cmp(1, 1 + 1e-6L) == -1, "cmp:容差外取符号");

    // ===== 2. p2 构造 / == / < / + - * / 与手算对照 =====
    p2 a = P(3, 4), b = P(-1, 2);
    CHECK(a.x == 3 && a.y == 4, "p2 聚合构造 {x, y}");
    CHECK(a == P(3, 4) && !(a == b), "operator== 相等/不等");
    CHECK(a == P(3 + eps / 4, 4) && a == P(3, 4 - eps / 4), "operator== 容差内视为相等");
    CHECK(!(a == P(3 + eps * 10, 4)), "operator== 容差外不等");
    CHECK(P(1, 2) < P(2, 0) && P(-1, 5) < P(0, -100), "operator< 先比 x");
    CHECK(P(1, 2) < P(1, 3) && !(P(1, 3) < P(1, 2)) && !(P(1, 2) < P(1, 2)), "operator< x 相等时比 y");
    CHECK(eqp(a + b, P(2, 6)) && eqp(b + a, P(2, 6)), "operator+ 与手算一致");
    CHECK(eqp(a - b, P(4, 2)) && eqp(b - a, P(-4, -2)), "operator- 与手算一致");
    CHECK(eqd(a * b, 5) && eqd(a * b, 3.0L * -1 + 4.0L * 2), "operator*(p2,p2) 点积(手算 3*(-1)+4*2 = 5)");
    CHECK(eqp(a * 2.0L, P(6, 8)) && eqp(2.0L * a, P(6, 8)) && eqp(a * -1.0L, P(-3, -4)), "operator* 数乘(db 在两侧)");
    CHECK(eqp(a / 2.0L, P(1.5L, 2.0L)) && eqp(a / -2.0L, P(-1.5L, -2.0L)), "operator/(p2,db)");
    CHECK(eqd(a.det(b), 10) && eqd(a.det(b), 3.0L * 2 - 4.0L * -1), "p2::det 叉积与手算一致(3*2-4*(-1) = 10)");
    CHECK(eqd(P(1, 0).det(P(0, 1)), 1) && eqd(P(0, 1).det(P(1, 0)), -1) && eqd(a.det(a), 0), "det 反对称/共线为 0");
    CHECK(eqd(P(1, 1).alpha(), pi / 4) && eqd(P(-1, 0).alpha(), pi) && eqd(P(0, -1).alpha(), -pi / 2), "alpha = atan2");
    CHECK(eqd(dis2(P(3, 4)), 25) && eqd(dis(P(3, 4)), 5) && eqd(dis(P(0, 0)), 0), "dis2/dis 长度");
    CHECK(eqd(dis2(P(-3, -4)), 25) && eqd(dis(P(3, 4) - P(3, 4)), 0), "dis2/dis 负数与零向量");
    CHECK(eqp(unit(P(3, 4)), P(0.6L, 0.8L)) && eqd(dis(unit(P(-3, 4))), 1), "unit 单位化");

    // ===== 3. r90 / rot =====
    CHECK(eqp(r90(P(3, 4)), P(-4, 3)) && eqp(r90(P(1, 0)), P(0, 1)), "r90 逆时针 90 度");
    CHECK(eqd(P(3, 4) * r90(P(3, 4)), 0) && eqd(dis2(r90(P(3, 4))), 25), "r90 正交且保长");
    CHECK(eqp(rot(P(1, 0), pi / 2), P(0, 1)) && eqp(rot(P(1, 0), pi), P(-1, 0)) && eqp(rot(P(0, 1), pi / 2), P(-1, 0)), "rot 特殊角");
    {
        int bad = 0;
        For(t, 1, 3000) {
            p2 v = P((db) rnd(-1000, 1000), (db) rnd(-1000, 1000));
            db ang = (db) rnd(-1000000, 1000000) / 100000;
            p2 r = rot(v, ang);
            if(!eqd(r.x, v.x * cos(ang) - v.y * sin(ang), 1e-9 * (1 + Abs(v.x)))) ++bad;
            if(!eqd(r.y, v.x * sin(ang) + v.y * cos(ang), 1e-9 * (1 + Abs(v.y)))) ++bad;
            if(!eqd(dis(r), dis(v), 1e-9 * (1 + dis(v)))) ++bad;
            if(!eqd(remainderl(r.alpha() - v.alpha() - ang, 2 * pi), 0, 1e-9)) ++bad;
        }
        CHECK(bad == 0, "rot 与手算公式、保长、角度加成一致(3000 组随机)");
    }

    // ===== 4. cross / crossop / inc(左转、右转、共线) =====
    seg L = {P(0, 0), P(2, 0)};
    CHECK(eqd(cross(L, P(1, 1)), 2) && crossop(L, P(1, 1)) == 1, "cross 左转为正");
    CHECK(eqd(cross(L, P(1, -1)), -2) && crossop(L, P(1, -1)) == -1, "cross 右转为负");
    CHECK(crossop(L, P(1, 0)) == 0 && crossop(L, P(5, 0)) == 0 && crossop(L, P(1, eps / 2)) == 0, "cross 共线为 0");
    CHECK(eqd(cross(L, P(1, 1)), L.dir().det(P(1, 1) - L.x)), "cross == dir().det(q - x)");
    CHECK(inc(L, P(9, 0)) && inc(L, P(0, 0)) && !inc(L, P(0, 1)), "inc 只看直线不看线段范围");
    CHECK(crossop({P(0, 0), P(1, 1)}, P(2, 2)) == 0 && crossop({P(0, 0), P(1, 1)}, P(2, 3)) == 1, "crossop 斜线段左转/共线");

    // ===== 5. isMid / ons / ons_s =====
    CHECK(isMid(0.0L, 0.5L, 1.0L) && isMid(1.0L, 0.5L, 0.0L), "isMid 区间内(两端顺序无关)");
    CHECK(isMid(0.0L, 0.0L, 1.0L) && isMid(0.0L, 1.0L, 1.0L) && isMid(2.0L, 2.0L, 2.0L), "isMid 端点算在区间内");
    CHECK(!isMid(0.0L, 1.5L, 1.0L) && !isMid(0.0L, -0.5L, 1.0L), "isMid 区间外为假");
    CHECK(isMid(P(0, 0), P(1, 1), P(2, 2)) && isMid(P(0, 0), P(1, 2), P(2, 2)), "isMid 点版:逐坐标落在区间内(即落在包围盒里)");
    CHECK(!isMid(P(0, 0), P(3, 2), P(2, 2)) && !isMid(P(0, 0), P(-1, 1), P(2, 2)), "isMid 点版:任一坐标出界即为假");
    CHECK(ons(L, P(1, 0)) && ons(L, P(0, 0)) && ons(L, P(2, 0)), "ons 线段内部与端点");
    CHECK(!ons(L, P(3, 0)) && !ons(L, P(-1, 0)) && !ons(L, P(1, eps * 100)), "ons 延长线与线外为假");
    CHECK(ons_s(L, P(1, 0)) && !ons_s(L, P(0, 0)) && !ons_s(L, P(2, 0)) && !ons_s(L, P(1, eps * 100)), "ons_s 严格在线段内部");
    CHECK(ons({P(0, 0), P(0, 3)}, P(0, 1.5)) && ons({P(1, 1), P(3, 3)}, P(2, 2)), "ons 竖直线段/斜线段");

    // ===== 6. proj / reflect / nearest(返回平方距离) =====
    seg X = {P(0, 0), P(4, 0)};
    CHECK(eqp(proj(X, P(1, 3)), P(1, 0)) && eqp(proj(X, P(-7, 3)), P(-7, 0)), "proj 垂足(含线段外)");
    CHECK(eqp(reflect(X, P(1, 3)), P(1, -3)) && eqp(reflect(X, P(5, -2)), P(5, 2)), "reflect 镜像");
    CHECK(eqd(nearest(X, P(1, 3)), 9) && eqd(nearest(X, P(-3, 1)), 10) && eqd(nearest(X, P(5, 1)), 2), "nearest 垂足在内/线段外(返回平方距离)");
    CHECK(eqd(nearest(X, P(0, 0)), 0) && eqd(nearest(X, P(4, 0)), 0), "nearest 端点距离为 0");
    {
        int bad = 0;
        For(t, 1, 30000) {
            p2 p1 = P((db) rnd(-1000, 1000), (db) rnd(-1000, 1000)), p2r = P((db) rnd(-1000, 1000), (db) rnd(-1000, 1000));
            if(dd(p1, p2r) < 4) continue;
            p2 q = P((db) rnd(-2000, 2000), (db) rnd(-2000, 2000));
            seg s = {p1, p2r};
            db ref = refNearest(p1, p2r, q);
            if(!eqd(nearest(s, q), ref, 1e-9 * (1 + ref))) ++bad;
            p2 h = proj(s, q);
            if(!eqd(dd(h, p1 + s.dir() * ((s.dir() * (q - p1)) / dis2(s.dir()))), 0, 1e-9 * (1 + dd(h, p1)))) ++bad;
            if(!eqd((q - h) * s.dir(), 0, 1e-6 * dis(s.dir()) * (1 + dis(q - p1)))) ++bad;
            if(!eqd(dd(reflect(s, q), h * 2 - q), 0, 1e-9)) ++bad;
            if(ons(s, q) != (ref < 1e-18 * (1 + dis2(s.dir())))) ++bad;
        }
        CHECK(bad == 0, "nearest/proj/reflect/ons 与独立参考一致(3 万组随机)");
    }

    // ===== 7. 线段相交 / 直线交点 =====
    CHECK(chkss({P(0, 0), P(2, 2)}, {P(0, 2), P(2, 0)}) && chkss_s({P(0, 0), P(2, 2)}, {P(0, 2), P(2, 0)}), "chkss/chkss_s 十字相交");
    CHECK(chkss({P(0, 0), P(1, 1)}, {P(1, 1), P(2, 0)}), "chkss 端点相接算相交");
    CHECK(!chkss_s({P(0, 0), P(1, 1)}, {P(1, 1), P(2, 0)}), "chkss_s 端点相接不算规范相交");
    CHECK(chkss({P(0, 0), P(2, 0)}, {P(1, 0), P(3, 0)}) && !chkss({P(0, 0), P(1, 0)}, {P(2, 0), P(3, 0)}), "chkss 共线重叠/共线分离");
    CHECK(!chkss({P(0, 0), P(1, 0)}, {P(0, 1), P(1, 1)}) && !chkss({P(0, 0), P(1, 1)}, {P(2, 2), P(3, 3)}), "chkss 平行/共线不重叠");
    seg M = {P(0, 1), P(2, 1)}, N = {P(0, 0), P(0, 2)}, L2 = {P(-5, 0), P(5, 0)};
    CHECK(!chkll(L, M) && !chkll(L, L) && chkll(L, N), "chkll 平行(含重合)为 0、相交为 1");
    CHECK(eqll(L, L2) && !eqll(L, M) && !eqll(L, N), "eqll 同一条直线/平行/相交");
    CHECK(eqp(isll(L, N), P(0, 0)), "isll 交点手算 x 轴 ∩ y 轴");
    CHECK(eqp(isll({P(0, 0), P(2, 2)}, {P(0, 2), P(2, 0)}), P(1, 1)), "isll 交点手算 两条对角线");
    CHECK(eqp(isll({P(1, 1), P(3, 1)}, {P(2, -1), P(2, 5)}), P(2, 1)), "isll 交点手算 水平 ∩ 垂直");
    {
        int bad = 0;
        For(t, 1, 30000) {
            p2 a1 = P((db) rnd(-10, 10), (db) rnd(-10, 10)), a2 = P((db) rnd(-10, 10), (db) rnd(-10, 10));
            p2 b1 = P((db) rnd(-10, 10), (db) rnd(-10, 10)), b2 = P((db) rnd(-10, 10), (db) rnd(-10, 10));
            if(dd(a1, a2) < 1 || dd(b1, b2) < 1) continue;
            if(chkss(seg{a1, a2}, seg{b1, b2}) != refChkss(a1, a2, b1, b2)) ++bad;
            int d1 = sign((a2 - a1).det(b1 - a1)), d2 = sign((a2 - a1).det(b2 - a1));
            int d3 = sign((b2 - b1).det(a1 - b1)), d4 = sign((b2 - b1).det(a2 - b1));
            if(chkss_s(seg{a1, a2}, seg{b1, b2}) != (d1 * d2 < 0 && d3 * d4 < 0)) ++bad;
            db den = (a2 - a1).det(b2 - b1);
            if(chkll(seg{a1, a2}, seg{b1, b2}) != (Abs(den) > 1e-9)) ++bad;
            if(Abs(den) > 1e-9) {
                p2 p = isll(seg{a1, a2}, seg{b1, b2});
                if(!eqd((a2 - a1).det(p - a1), 0, 1e-6 * (1 + dis(a2 - a1)))) ++bad;
                if(!eqd((b2 - b1).det(p - b1), 0, 1e-6 * (1 + dis(b2 - b1)))) ++bad;
            } else if(eqll(seg{a1, a2}, seg{b1, b2}) != (Abs((b1 - a1).det(a2 - a1)) < 1e-9)) ++bad;
        }
        CHECK(bad == 0, "chkss/chkss_s/chkll/isll/eqll 与独立参考一致(3 万组随机)");
    }

    // ===== 8. disss(返回平方距离) =====
    CHECK(eqd(disss({P(0, 0), P(1, 0)}, {P(0, 1), P(1, 1)}), 1), "disss 平行线段间隔为 1");
    CHECK(eqd(disss({P(0, 0), P(2, 2)}, {P(0, 2), P(2, 0)}), 0), "disss 相交为 0");
    CHECK(eqd(disss({P(0, 0), P(1, 0)}, {P(3, 2), P(4, 2)}), 8), "disss 端点间最短(2^2+2^2 = 8)");
    CHECK(eqd(disss({P(0, 0), P(1, 0)}, {P(3, 4), P(4, 4)}), 20), "disss 端点间最短(2^2+4^2 = 20)");
    {
        int bad = 0;
        For(t, 1, 30000) {
            p2 a1 = P((db) rnd(-20, 20), (db) rnd(-20, 20)), a2 = P((db) rnd(-20, 20), (db) rnd(-20, 20));
            p2 b1 = P((db) rnd(-20, 20), (db) rnd(-20, 20)), b2 = P((db) rnd(-20, 20), (db) rnd(-20, 20));
            if(dd(a1, a2) < 1 || dd(b1, b2) < 1) continue;
            db ref = refDisss(a1, a2, b1, b2);
            if(!eqd(disss(seg{a1, a2}, seg{b1, b2}), ref, 1e-9 * (1 + ref))) ++bad;
        }
        CHECK(bad == 0, "disss 与独立参考一致(3 万组随机)");
    }

    // ===== 9. area / contain =====
    p2 sq[4] = {P(0, 0), P(2, 0), P(2, 2), P(0, 2)};
    p2 sqr[4] = {P(0, 0), P(0, 2), P(2, 2), P(2, 0)};
    CHECK(eqd(area(4, sq), 4) && eqd(area(4, sqr), -4), "area 带号面积(逆时针为正、顺时针为负)");
    {
        int bad = 0;
        For(x, -2, 6) For(y, -2, 6) {
            int got = contain(4, sq, P((db) x / 2, (db) y / 2));
            bool in = x >= 0 && x <= 4 && y >= 0 && y <= 4;
            bool bd = in && (x == 0 || x == 4 || y == 0 || y == 4);
            int want = bd ? 1 : (in ? 2 : 0);
            if(got != want) ++bad;
        }
        CHECK(bad == 0, "contain 正方形网格 0(外)/1(边界)/2(内)全部正确");
    }
    {
        int bad = 0;
        For(t, 1, 200) {
            vector<p2> v((int) rnd(3, 10));
            ForD(i, 0, v.size()) v[i] = P((db) rnd(-40, 40), (db) rnd(-40, 40));
            p2 hb[32];
            int k = convex_hull(v.size(), v.data(), hb);
            if(k < 3) continue;
            vector<p2> h(hb, hb + k);
            if(!eqd(Abs(area(k, h.data())), polyArea(h))) ++bad;
            For(it, 1, 200) {
                p2 q = P((db) rnd(-60, 60), (db) rnd(-60, 60));
                if(contain(k, h.data(), q) != refContain(h, q)) ++bad;
            }
        }
        CHECK(bad == 0, "area/contain 与独立参考一致(200 个随机凸多边形 × 200 点)");
    }

    // ===== 10. convex_hull / convex_diameter / convex_cut =====
    {
        int badh = 0, badd = 0;
        For(t, 1, 3000) {
            int n = (int) rnd(2, 14);
            vector<p2> v;
            while((int) v.size() < n) {  // 造互不相同的点(本模板与上游板子都按“去重点”使用)
                p2 q = P((db) rnd(-50, 50), (db) rnd(-50, 50));
                bool dup = false;
                ForD(i, 0, v.size()) if(eqp(v[i], q)) dup = true;
                if(!dup) v.push_back(q);
            }
            p2 a2[64], hb[64];
            ForD(i, 0, n) a2[i] = v[i];
            int k = convex_hull(n, a2, hb);
            vector<p2> mine(hb, hb + k), ref = refHull(v);
            if(mine.size() != ref.size()) { ++badh; badd++; continue; }
            int off = -1;
            ForD(i, 0, mine.size()) if(eqp(mine[i], ref[0])) off = i;
            if(off < 0) { ++badh; badd++; continue; }
            bool same = true;
            ForD(i, 0, ref.size()) if(!eqp(mine[(off + i) % mine.size()], ref[i])) same = false;
            if(!same) ++badh;
            ForD(i, 0, mine.size()) {  // 逆时针凸性
                db c = cross(seg{mine[i], mine[(i + 1) % mine.size()]}, mine[(i + 2) % mine.size()]);
                if(c < -1e-6) ++badh;
            }
            ForD(i, 0, n) ForD(j, 0, mine.size())  // 所有输入点都在凸包内(含边界)
                if(cross(seg{mine[j], mine[(j + 1) % mine.size()]}, v[i]) < -1e-6) ++badh;
            db bf = 0;
            ForD(i, 0, mine.size()) ForD(j, 0, mine.size()) bf = max(bf, dd(mine[i], mine[j]));
            if(!eqd(convex_diameter(mine.size(), mine.data()), bf, 1e-9 * (1 + bf))) ++badd;
        }
        CHECK(badh == 0, "convex_hull 与独立参考一致(3000 组随机):顶点、逆时针顺序、凸性、包含所有点");
        CHECK(badd == 0, "convex_diameter 与暴搜最远点对一致(返回平方距离,3000 组)");
    }
    {
        // nos=0 去掉边上的共线点,nos=1 保留
        vector<p2> pts = {P(0, 0), P(1, 0), P(2, 0), P(3, 0), P(3, 3), P(0, 3), P(0, 1), P(1, 1)};
        vector<p2> t0 = pts, t1 = pts;
        p2 h0[32], h1[32];
        int k0 = convex_hull(t0.size(), t0.data(), h0, 0), k1 = convex_hull(t1.size(), t1.data(), h1, 1);
        CHECK(k0 == 4 && eqp(h0[0], P(0, 0)) && eqp(h0[1], P(3, 0)) && eqp(h0[2], P(3, 3)) && eqp(h0[3], P(0, 3)), "convex_hull nos=0 去掉共线点(只剩 4 个角点)");
        CHECK(k1 == 7, "convex_hull nos=1 保留共线点(4 角点 + 3 个边上点)");
    }
    {
        int bad = 0;
        For(t, 1, 1000) {
            vector<p2> v((int) rnd(3, 10));
            ForD(i, 0, v.size()) v[i] = P((db) rnd(-30, 30), (db) rnd(-30, 30));
            p2 hb[32];
            int k = convex_hull(v.size(), v.data(), hb);
            if(k < 3) continue;
            vector<p2> h(hb, hb + k);
            seg q = {P((db) rnd(-60, 60), (db) rnd(-60, 60)), P((db) rnd(-60, 60), (db) rnd(-60, 60))};
            if(dd(q.x, q.y) < 1) continue;
            vector<p2> got = convex_cut(k, h.data(), q), want = refCut(h, q);
            if(got.size() != want.size()) { ++bad; continue; }
            ForD(i, 0, want.size()) if(!eqp(got[i], want[i], 1e-6)) ++bad;
            if(!eqd(polyArea(got), polyArea(want), 1e-6 * (1 + polyArea(want)))) ++bad;
        }
        CHECK(bad == 0, "convex_cut 与独立参考半平面裁剪一致(1000 组随机)");
    }

    // ===== 11. circumcircle_diameter / circumcenter =====
    {
        auto c = circumcenter(P(0, 0), P(2, 0), P(1, 1));
        CHECK(eqd(circumcircle_diameter(P(0, 0), P(2, 0), P(1, 1)), 2) && eqp(P(c[0], c[1]), P(1, 0)), "外接圆:直角三角形斜边为直径、圆心在斜边中点");
    }
    CHECK(circumcircle_diameter(P(0, 0), P(1, 0), P(2, 0)) == -1, "共线三点返回 -1");
    {
        int bad = 0;
        For(t, 1, 30000) {
            p2 A = P((db) rnd(-50, 50), (db) rnd(-50, 50)), B = P((db) rnd(-50, 50), (db) rnd(-50, 50)), C = P((db) rnd(-50, 50), (db) rnd(-50, 50));
            db ar2 = Abs((B - A).det(C - A));
            if(ar2 < 1) continue;
            db R = dis(A - B) * dis(B - C) * dis(C - A) / (2 * ar2);  // abc/(4Δ),4Δ = 2|det|
            if(!eqd(circumcircle_diameter(A, B, C), 2 * R, 1e-9 * (1 + R))) ++bad;
            auto cc = circumcenter(A, B, C);
            p2 O = P(cc[0], cc[1]);
            if(!eqd(dis(O - A), R, 1e-9 * (1 + R)) || !eqd(dis(O - B), R, 1e-9 * (1 + R)) || !eqd(dis(O - C), R, 1e-9 * (1 + R))) ++bad;
        }
        CHECK(bad == 0, "circumcircle_diameter(= 2R)与 circumcenter 到三点等距(3 万组随机三角形)");
    }

    // ===== 12. 边界观察(不计入失败):n <= 1 时 convex_hull 不写 b =====
    {
        p2 only[1] = {P(7, 8)};
        p2 out[4] = {P(-1, -1), P(-1, -1), P(-1, -1), P(-1, -1)};
        if(convex_hull(1, only, out) == 1 && !eqp(out[0], only[0]))
            printf("  [note] convex_hull(n=1) 返回 1 但不往 b 里写点(b[0] 仍是哨兵 (-1,-1)),n<=1 需调用方自己处理\n");
    }

    // ===== 13. operator+= (geo.cpp:32,正确) =====
    {
        int bad = 0;
        For(t, 1, 20000) {
            p2 x = P((db) rnd(-1000, 1000), (db) rnd(-1000, 1000)), y = P((db) rnd(-1000, 1000), (db) rnd(-1000, 1000));
            p2 x0 = x, r = (x += y);
            if(!eqp(x, x0 + y) || !eqp(r, x0 + y)) ++bad;
        }
        CHECK(bad == 0, "operator+=(p2&,p2) 语义正确:a += b 后 a == a + b(geo.cpp:32,2 万组随机)");
    }

    // ===== 14. 复合赋值 -= 与 /=:当前 geo.cpp:33-34 把这两个都写成了加法 =====
    // geo.cpp:33  p2 operator-=(p2 &x, p2 y) { return x = x + y; }   应为 x = x - y
    // geo.cpp:34  p2 operator/=(p2 &x, p2 y) { return x = x + y; }   应为 x = x / y(第二参数应为 db)
    bool minus_ok = false, diveq_ok = false;
    {  // --- operator-= ---
        p2 x = P(3, 4), y = P(1, 1), x0 = x, want = x0 - y;
        p2 r = (x -= y);
        minus_ok = eqp(x, want) && eqp(r, want);
        int same_add = 0, same_sub = 0;
        For(t, 1, 200) {
            p2 u = P((db) rnd(-100, 100), (db) rnd(-100, 100)), v = P((db) rnd(-100, 100), (db) rnd(-100, 100));
            p2 w = u;
            w -= v;
            if(eqp(w, u + v)) ++same_add;
            if(eqp(w, u - v)) ++same_sub;
        }
        printf("  %s geo.cpp:33 operator-=(p2 &x, p2 y)%s\n", minus_ok ? "[ok] " : "[BUG]",
               minus_ok ? "" : " 当前源码是 x = x + y(加法),应为 x = x - y");
        printf("        ");
        pr(x0);
        printf(" -= ");
        pr(y);
        printf("  ->  实测 a = ");
        pr(x);
        printf(",期望 a = ");
        pr(want);
        printf(" (= a - b);返回值 ");
        pr(r);
        printf(" 与实测 a 一致\n");
        printf("        200 组随机用例:%d 组结果等于 a + b(期望 0 组),%d 组等于 a - b(期望 200 组)\n", same_add, same_sub);
    }
    {  // --- operator/= ---
        const bool has_db = has_diveq_db<p2>::value, has_p2 = has_diveq_p2<p2>::value;
        p2 x = P(6, 9), y = P(2, 2), out, out2;
        bool ran_db = run_diveq_db(x, 2.0L, out);
        bool db_ok = ran_db && eqp(out, x / 2.0L);
        bool ran_p2 = run_diveq_p2(x, y, out2);
        int same_add = 0;
        if(ran_p2) {  // 只能通过模板助手调用:修好后 (p2&, p2) 重载不存在,直接写 w /= v 会编译失败
            For(t, 1, 200) {
                p2 u = P((db) rnd(-100, 100), (db) rnd(-100, 100)), v = P((db) rnd(-100, 100), (db) rnd(-100, 100)), w;
                run_diveq_p2(u, v, w);
                if(eqp(w, u + v)) ++same_add;
            }
        }
        diveq_ok = has_db && db_ok;
        printf("  %s geo.cpp:34 operator/=(p2 &x, ?)%s\n", diveq_ok ? "[ok] " : "[BUG]",
               diveq_ok ? "" : " 当前源码是 (p2&, p2) 版且写成了 x = x + y(加法)");
        printf("        重载探测:(p2&, db) 版存在 = %d(期望 1),(p2&, p2) 版存在 = %d(期望 0)\n", (int) has_db, (int) has_p2);
        if(ran_db) {
            printf("        a /= 2.0 实测:");
            pr(x);
            printf("  ->  ");
            pr(out);
            printf(",期望 ");
            pr(x / 2.0L);
            printf(" (= a / 2)\n");
        } else {
            printf("        (p2&, db) 重载缺失:`a /= 2.0` 根本编译不过 —— error: no match for 'operator/='\n");
        }
        if(ran_p2) {
            printf("        a /= b 实测:");
            pr(x);
            printf(" /= ");
            pr(y);
            printf("  ->  ");
            pr(out2);
            printf(",其中 %d/200 组随机用例等于 a + b(加法,期望 0 组)\n", same_add);
        }
    }
    if(minus_ok && diveq_ok) {
        printf("  ==> 结论:+= / -= /= 三个复合赋值运算符语义都正确\n");
    } else {
        printf("  ==> 结论:+= 正确,但 -= 与 /= 复现了已知 bug(geo.cpp:33-34 都写成 x = x + y),本 check 因此 FAIL\n");
        printf("      (其余接口在上面全部 [ok],只差这两个运算符)\n");
    }
    CHECK(minus_ok && diveq_ok, "已知 bug:geo.cpp:33 operator-= 与 geo.cpp:34 operator/= 实现成了加法(应分别为 x - y 与 x / y)");

    PASSED("geo");
}
