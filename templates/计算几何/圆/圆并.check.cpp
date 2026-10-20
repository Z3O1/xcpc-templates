// 圆并 自测:自适应积分/分段周长参考、500组随机，旋转换根/隐藏根、退化/洞/跨角/极端缩放/远簇/2000圆
#include "../../_check_base.hpp"
#include "../基础/二维向量.cpp"
#include "../基础/直线与线段.cpp"
#include "圆基础.cpp"
#include "圆并.cpp"

// 失败打印输入圆；相对误差按单位化后的面积、长度度量。
void need(bool good, const vector<c2> &a, const char *what) {
    if(good) return;
    printf("FAIL %s n=%zu\n", what, a.size());
    for(auto c : a) printf("(%Lg,%Lg;r=%Lg) ", c.o.x, c.o.y, c.r);
    puts("");
    exit(1);
}
// 固定x处，各圆纵向区间合并后的总长；不使用角扫描或生产圆求交。
db height(const vector<c2> &a, db x) {
    vector<array<db, 2>> v;
    for(auto c : a) {
        db d = abs(x - c.o.x);
        if(d >= c.r) continue;
        db h = sqrt((c.r - d) * (c.r + d));
        v.push_back({c.o.y - h, c.o.y + h});
    }
    sort(v.begin(), v.end());
    db sum = 0, end = -numeric_limits<db>::infinity();
    for(auto p : v) {
        sum += max(db(0), p[1] - max(p[0], end));
        end = max(end, p[1]);
    }
    return sum;
}
// 自适应Simpson参考，限制递归深度，积分函数与被测算法独立。
db integrate(const function<db(db)> &f, db l, db r, db fl, db fm, db fr, db s, db tol, int depth) {
    db m = l + (r - l) / 2, x = l + (m - l) / 2, y = m + (r - m) / 2;
    db fx = f(x), fy = f(y);
    db left = (m - l) * (fl + 4 * fx + fm) / 6, right = (r - m) * (fm + 4 * fy + fr) / 6;
    if(depth == 0 || abs(left + right - s) < 15 * tol) return left + right + (left + right - s) / 15;
    return integrate(f, l, m, fl, fx, fm, left, tol / 2, depth - 1) +
           integrate(f, m, r, fm, fy, fr, right, tol / 2, depth - 1);
}
// 用所有圆的左右端点分段，纵切积分得到并面积。
db reference_area(const vector<c2> &a) {
    vector<db> xs;
    for(auto c : a) if(c.r > 0) { xs.push_back(c.o.x - c.r); xs.push_back(c.o.x + c.r); }
    sort(xs.begin(), xs.end());
    xs.erase(unique(xs.begin(), xs.end()), xs.end());
    function<db(db)> f = [&](db x) { return height(a, x); };
    db answer = 0;
    for(int i = 1; i < (int) xs.size(); ++i) {
        db l = xs[i - 1], r = xs[i], fl = f(l), fm = f((l + r) / 2), fr = f(r);
        answer += integrate(f, l, r, fl, fm, fr, (r - l) * (fl + 4 * fm + fr) / 6, 1e-11L, 24);
    }
    return answer;
}
// 独立解析求边界交点，再逐角段取中点判断暴露；不用生产半角公式或覆盖区间合并。
db reference_perimeter(const vector<c2> &a) {
    db answer = 0;
    for(int i = 0; i < (int) a.size(); ++i) {
        auto c = a[i];
        if(c.r == 0) continue;
        bool hidden = false;
        vector<db> angles{0, 2 * pi};
        for(int j = 0; j < (int) a.size(); ++j) if(j != i) {
            auto b = a[j];
            p2 v = b.o - c.o;
            db d = hypot(v.x, v.y);
            if(d + c.r <= b.r && (c.r < b.r || j < i)) { hidden = true; break; }
            if(d == 0 || d >= c.r + b.r || d <= abs(c.r - b.r)) continue;
            db x = (c.r * c.r - b.r * b.r + d * d) / (2 * d);
            db y = sqrt(max(db(0), c.r * c.r - x * x));
            p2 u{v.x / d, v.y / d};
            for(int sign : {-1, 1}) {
                p2 q{u.x * x - u.y * y * sign, u.y * x + u.x * y * sign};
                db theta = atan2(q.y, q.x);
                if(theta < 0) theta += 2 * pi;
                angles.push_back(theta);
            }
        }
        if(hidden) continue;
        sort(angles.begin(), angles.end());
        for(int k = 1; k < (int) angles.size(); ++k) {
            db m = (angles[k - 1] + angles[k]) / 2;
            p2 offset{c.r * cos(m), c.r * sin(m)};
            bool outside = true;
            for(int j = 0; j < (int) a.size(); ++j) if(j != i && a[j].r > 0) {
                if(a[j].r == c.r && a[j].o.x == c.o.x && a[j].o.y == c.o.y) continue;
                p2 v = c.o - a[j].o + offset;
                if(v.x * v.x + v.y * v.y < a[j].r * a[j].r) outside = false;
            }
            if(outside) answer += c.r * (angles[k] - angles[k - 1]);
        }
    }
    return answer;
}
// 与面积数值积分、独立周长参考及基本界对照；兼测输入不变。
void verify(vector<c2> a) {
    auto before = a, normalized = a;
    auto got = cir_union(a);
    db scale = 0;
    for(auto c : a) scale = max(scale, c.r);
    if(scale == 0) { need(got[0] == 0 && got[1] == 0, a, "空集/全点圆"); return; }
    for(auto &c : normalized) c.o = (c.o - a[0].o) / scale, c.r /= scale;
    db s = reference_area(normalized), p = reference_perimeter(normalized);
    db gs = got[0] / scale / scale, gp = got[1] / scale;
    if(abs(gs - s) > 3e-7L * (1 + s) || abs(gp - p) > 2e-10L * (1 + p))
        printf("got=(%Lg,%Lg) want=(%Lg,%Lg)\n", gs, gp, s, p);
    need(abs(gs - s) < 3e-7L * (1 + s) && abs(gp - p) < 2e-10L * (1 + p), a, "独立面积/周长参考");
    db sum = 0, maxarea = 0, length = 0;
    for(auto c : normalized) { sum += pi * c.r * c.r; maxarea = max(maxarea, pi * c.r * c.r); length += 2 * pi * c.r; }
    need(gs >= maxarea - 1e-10L && gs <= sum + 1e-10L && gp >= 0 && gp <= length + 1e-10L &&
         gp * gp >= 4 * pi * gs - 1e-8L, a, "面积/周长界及等周性质");
    for(int i = 0; i < (int) a.size(); ++i)
        need(a[i].o.x == before[i].o.x && a[i].o.y == before[i].o.y && a[i].r == before[i].r, a, "输入不变");
}
// 精确或解析结果，按尺度归一化后比较，包含不能用积分参考的大平移/悬殊半径情形。
void analytic(const vector<c2> &a, db s, db p, db scale = 1) {
    auto got = cir_union(a);
    need(abs(got[0] / scale / scale - s) < 2e-10L * (1 + abs(s)) &&
         abs(got[1] / scale - p) < 2e-10L * (1 + abs(p)), a, "解析回归");
}
// 随机数值对拍与多种退化、尺度、拓扑和大规模解析回归。
int main() {
    mt19937_64 gen(2026101603);
    verify({}); verify({{{0, 0}, 0}});
    for(int t = 0; t < 500; ++t) {
        vector<c2> a(gen() % 7);
        for(auto &c : a) c = {{db(int(gen() % 81) - 40) / 10, db(int(gen() % 81) - 40) / 10}, db(gen() % 30 + 1) / 10};
        if(t % 5 == 0 && !a.empty()) a.push_back(a[0]);
        if(t % 7 == 0) a.push_back({{0, 0}, 0});
        verify(a);
        if(t % 10 == 0 && !a.empty()) {
            auto expected = cir_union(a);
            shuffle(a.begin(), a.end(), gen);
            auto got = cir_union(a);
            need(abs(got[0] - expected[0]) < 1e-10L && abs(got[1] - expected[1]) < 1e-10L, a, "置换不变");
            db scale = powl(10, t % 25 - 12);
            for(auto &c : a) c.o = rot(c.o, 0.731L) * scale + p2{17 * scale, -23 * scale}, c.r *= scale;
            analytic(a, expected[0], expected[1], scale);
        }
    }
    db two_area = 4 * pi / 3 + sqrt(db(3)) / 2, two_perimeter = 8 * pi / 3;
    // 同一解析并集旋转、平移并换根；覆盖小圆也可能成为连通块原点。
    for(db scale : {db(1e-18L), db(1e-12L), db(1), db(1e12L), db(1e40L)})
        for(int angle = 0; angle < 16; ++angle) for(bool hidden_root : {false, true}) {
            p2 origin{16 * scale, -16 * scale}, v = rot({scale, 0}, angle * pi / 8);
            vector<c2> cluster{{origin, scale}, {origin + v, scale}};
            if(hidden_root) cluster.push_back({origin + v / 2, scale / 4});
            for(int i = 0; i < (int) cluster.size(); ++i) {
                analytic(cluster, two_area, two_perimeter, scale);
                rotate(cluster.begin(), cluster.begin() + 1, cluster.end());
            }
        }
    for(db scale : {db(1e-18L), db(1e-12L), db(1), db(1e12L), db(1e40L)})
        analytic({{{0, 0}, scale}, {{scale, 0}, scale}}, two_area, two_perimeter, scale);
    analytic({{{0, 0}, 1}, {{0, 0}, 1}, {{0, 0}, 0}}, pi, 2 * pi);
    analytic({{{0, 0}, 3}, {{1, 0}, 1}, {{2, 0}, 1}}, 9 * pi, 6 * pi);
    analytic({{{0, 0}, 1}, {{2, 0}, 1}}, 2 * pi, 4 * pi);
    analytic({{{0, 0}, 1}, {{0, 0}, 0}, {{10, 10}, 0}}, pi, 2 * pi);
    analytic({{{0, 0}, 1}, {{1, 0}, 1}, {{0.5L, sqrt(db(3)) / 2}, 1}}, 3 * pi / 2 + sqrt(db(3)), 3 * pi);
    db side = 1.8L, lens = 2 * acos(side / 2) - side * sqrt(4 - side * side) / 2;
    analytic({{{0, 0}, 1}, {{side, 0}, 1}, {{side / 2, side * sqrt(db(3)) / 2}, 1}},
             3 * pi - 3 * lens, 6 * pi - 12 * acos(side / 2));  // 中间有洞，周长包含内边界。
    for(db T : {ldexp(db(1), 40), db(1e18L)})
        analytic({{{0, 0}, 1}, {{1, 0}, 1}, {{T, T}, 1}, {{T + 1, T}, 1}}, 2 * two_area, 2 * two_perimeter);
    // 小圆中心恰在大圆边界；不能被 d+r 的舍入误判为完全包含/外离。
    auto unequal = cir_union({{{0, 0}, 1e18L}, {{1e18L, 0}, 1}});
    need(abs(unequal[0] / 1e18L / 1e18L - pi) < 1e-14L &&
         abs(unequal[1] / 1e18L - 2 * pi) < 1e-14L, {}, "悬殊半径");
    vector<c2> a(2000);
    for(int i = 0; i < 2000; ++i) a[i] = {{3.L * i, 0}, 1};
    analytic(a, 2000 * pi, 4000 * pi);
    fill(a.begin(), a.end(), c2{{0, 0}, 1});
    analytic(a, pi, 2 * pi);
    a.resize(1200);
    for(int i = 0; i < (int) a.size(); ++i) a[i] = {{cos(2 * pi * i / a.size()), sin(2 * pi * i / a.size())}, 1};
    analytic(a, 2 * pi + a.size() * sin(2 * pi / a.size()), 4 * pi);
    PASSED("圆并");
}
