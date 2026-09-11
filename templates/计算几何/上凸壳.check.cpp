// 上凸壳 自测:upper_hull 与「翻转 y 后求标准下凸壳」的独立参考对拍 + 性质断言
//
// 参考实现(文件内自写,不用模板):
//   refUpper(v):把每个点按 y 取反,套一遍标准 Andrew 下凸链(只留严格左转),再把 y 翻回来 ——
//     得到的正是上半凸链;它与「排序 + 单调栈求右转链」是两套等价但写法不同的实现。
//   refHullUp(v):另一条独立路线 —— 先求完整逆时针凸包(Andrew),再从「最右上角的点」沿凸包
//     往「最左上角的点」走,取上半圈顶点;两条参考互相校验后再与模板比。
//   性质级断言:① 链上相邻三点全右转(叉积 <= 0,严格时 < 0);② 所有输入点都在链上或链下方;
//     ③ 首末点是最左上的点与最右上的点;④ 链上每个点都来自输入;⑤ x 严格递增;⑥ 顶点数 <= n。
//
// 用例规模:
//   1 解析用例:正方形(上边界两个角)、矩形顶边塞满共线点(nos=0 vs nos=1)、全共线、
//     全同点、n=0/1/2/3、竖直边(同一 x 多点只留最高)、圆上点、凹点集
//   2 4000 组 n=1..40 小整数坐标(±12,含重复点/共线点)-> 两条参考逐位对照 + 全部性质断言
//   3 3000 组 n=1..60 大坐标(±1e6)只做参考对照与性质断言
//   4 3000 组「上边界共线点很多」的构造用例:nos=0/nos=1 的语义各自断言
//   5 1000 组 ±1e9 大坐标(整数,叉积精确)-> 参考逐位对照
#include "../_check_base.hpp"
#include "geo.cpp"
#include "上凸壳.cpp"

static p2 P(db x, db y) { return {x, y}; }
static db Abs(db x) { return x < 0 ? -x : x; }
static void pr(p2 a) { printf("(%g,%g)", (double) a.x, (double) a.y); }
static void prv(const char *s, const vector<p2> &v) {
    printf("      %s(%zu):", s, v.size());
    ForD(i, 0, (int) min<size_t>(v.size(), 30)) printf(" "), pr(v[i]);
    printf("%s\n", v.size() > 30 ? " ..." : "");
}
static vector<p2> runUpper(const vector<p2> &v, bool nos, int &k) {  // 跑模板并取出结果
    static p2 a[4096], b[4096];
    ForD(i, 0, (int) v.size()) a[i] = v[i];
    k = upper_hull(v.size(), a, b, nos);
    return vector<p2>(b, b + k);
}

// ---------- 独立参考实现 ----------
static vector<p2> topPerX(const vector<p2> &v) {  // 「每个 x 只留最高的点」是上凸壳的定义的一部分,两条参考都先做这一步
    vector<p2> w = v;
    sort(w.begin(), w.end(), [](const p2 &u, const p2 &t) { return u.x != t.x ? u.x < t.x : u.y > t.y; });
    vector<p2> r;
    ForD(i, 0, (int) w.size()) if(r.empty() || r.back().x != w[i].x) r.push_back(w[i]);
    return r;
}
static vector<p2> refLower(vector<p2> v) {  // 标准 Andrew 下凸链(只留严格左转)
    sort(v.begin(), v.end(), [](const p2 &u, const p2 &w) { return u.x != w.x ? u.x < w.x : u.y < w.y; });
    vector<p2> h;
    ForD(i, 0, (int) v.size()) {
        while(h.size() > 1 && sign((h[h.size() - 1] - h[h.size() - 2]).det(v[i] - h[h.size() - 2])) <= 0) h.pop_back();
        h.push_back(v[i]);
    }
    return h;
}
static vector<p2> refUpper(const vector<p2> &v) {  // y 取反 -> 标准下凸链 -> 翻回来 = 上半凸链
    vector<p2> w;
    for(p2 q : topPerX(v)) w.push_back(P(q.x, -q.y));
    vector<p2> h = refLower(w);
    ForD(i, 0, (int) h.size()) h[i].y = -h[i].y;
    return h;
}
static vector<p2> refHull(vector<p2> v) {  // 完整逆时针凸包(严格去共线)
    sort(v.begin(), v.end(), [](const p2 &u, const p2 &w) { return u.x != w.x ? u.x < w.x : u.y < w.y; });
    v.erase(unique(v.begin(), v.end(), [](p2 a, p2 b) { return a.x == b.x && a.y == b.y; }), v.end());
    int n = v.size();
    if(n <= 2) return v;
    vector<p2> h;
    ForD(i, 0, n) {
        while(h.size() > 1 && sign((h[h.size() - 1] - h[h.size() - 2]).det(v[i] - h[h.size() - 2])) <= 0) h.pop_back();
        h.push_back(v[i]);
    }
    int t = h.size();
    rFor(i, n - 2, 0) {
        while((int) h.size() > t && sign((h[h.size() - 1] - h[h.size() - 2]).det(v[i] - h[h.size() - 2])) <= 0) h.pop_back();
        h.push_back(v[i]);
    }
    h.pop_back();
    return h;
}
static vector<p2> refHullUp(const vector<p2> &v) {  // 从完整凸包上截出上半圈:从最左上点逆着走到最右上点
    vector<p2> h = refHull(topPerX(v));
    if(h.size() <= 1) return h;
    int k = h.size(), jl = 0, jr = 0;  // jl:最左上(min x 中 max y),jr:最右上(max x 中 max y)
    ForD(i, 1, k) {
        if(h[i].x < h[jl].x || (h[i].x == h[jl].x && h[i].y > h[jl].y)) jl = i;
        if(h[i].x > h[jr].x || (h[i].x == h[jr].x && h[i].y > h[jr].y)) jr = i;
    }
    vector<p2> up;
    if(k == 2) return h;
    for(int i = jl;; i = (i + k - 1) % k) {  // 逆时针倒着走 = 沿上半圈
        up.push_back(h[i]);
        if(i == jr) break;
        if((int) up.size() > k) break;
    }
    return up;
}
// 一条链上的性质断言(不依赖参考)
static int chainProps(const vector<p2> &v, const vector<p2> &h, bool nos, db scale) {
    if(h.empty()) {
        printf("      链为空但输入非空\n");
        return 0;
    }
    db tol = 1e-9L * scale;
    ForD(i, 0, (int) h.size()) {
        bool found = false;
        ForD(j, 0, (int) v.size()) if(h[i] == v[j]) found = true;
        if(!found) {
            printf("      链上第 %d 个点 不是输入点\n", i), pr(h[i]), printf("\n");
            return 0;
        }
    }
    ForD(i, 1, (int) h.size()) if(cmp(h[i - 1].x, h[i].x) >= 0) {
        printf("      链的 x 不是严格递增:第 %d 个点 ", i - 1), pr(h[i - 1]), printf(" -> "), pr(h[i]), printf("\n");
        return 0;
    }
    ForD(i, 0, (int) h.size() - 2) {
        db c = cross({h[i], h[i + 1]}, h[i + 2]);
        if(nos ? c > tol : c >= -tol) {
            printf("      第 %d 个点处不是右转:叉积 %.17Lg(nos=%d)\n", i + 1, c, (int) nos);
            prv("链", h);
            return 0;
        }
    }
    ForD(i, 0, (int) h.size() - 1) ForD(j, 0, (int) v.size()) {  // 所有输入点在链上或链下方
        if(cmp(v[j].x, h[i].x) < 0 || cmp(v[j].x, h[i + 1].x) > 0) continue;
        db c = cross({h[i], h[i + 1]}, v[j]);
        if(c > tol * (1 + dis(h[i + 1] - h[i]))) {
            printf("      输入点 "), pr(v[j]), printf(" 在链段 "), pr(h[i]), printf("-"), pr(h[i + 1]), printf(" 上方(叉积 %.17Lg)\n", c);
            return 0;
        }
    }
    db mnx = 1e30L, mxx = -1e30L;
    ForD(i, 0, (int) v.size()) mnx = min(mnx, v[i].x), mxx = max(mxx, v[i].x);
    db mxy0 = -1e30L, mxy1 = -1e30L;
    ForD(i, 0, (int) v.size()) {
        if(!cmp(v[i].x, mnx)) mxy0 = max(mxy0, v[i].y);
        if(!cmp(v[i].x, mxx)) mxy1 = max(mxy1, v[i].y);
    }
    if(!(h.front() == P(mnx, mxy0)) || !(h.back() == P(mxx, mxy1))) {
        printf("      首/末点不是最左上/最右上(最左上 ");
        pr(P(mnx, mxy0)), printf(",最右上 "), pr(P(mxx, mxy1)), printf(";实测首 "), pr(h.front()), printf(" 末 "), pr(h.back()), printf(")\n");
        return 0;
    }
    return 1;
}
static bool sameChain(const vector<p2> &x, const vector<p2> &y) {
    if(x.size() != y.size()) return false;
    ForD(i, 0, (int) x.size()) if(!(x[i] == y[i])) return false;
    return true;
}
static int oneCase(const vector<p2> &v, int verbose) {
    int k0 = 0, k1 = 0;
    vector<p2> h0 = runUpper(v, 0, k0), h1 = runUpper(v, 1, k1);
    vector<p2> r0 = refUpper(v), r1 = refHullUp(v);
    db mx = 1;
    ForD(i, 0, (int) v.size()) mx = max(mx, max(Abs(v[i].x), Abs(v[i].y)));
    if(!sameChain(h0, r0) || !sameChain(r0, r1)) {
        printf("      nos=0 的链与参考不一致(模板 %zu 点,翻转法 %zu 点,凸包截取 %zu 点)\n", h0.size(), r0.size(), r1.size());
        prv("输入", v), prv("模板", h0), prv("翻转法", r0), prv("凸包截取", r1);
        return 0;
    }
    if(!chainProps(v, h0, 0, mx) || !chainProps(v, h1, 1, mx)) {
        if(verbose) prv("输入", v);
        return 0;
    }
    if(h1.size() < h0.size()) {
        printf("      nos=1 的链比 nos=0 还短(%zu < %zu)\n", h1.size(), h0.size());
        return 0;
    }
    return 1;
}

int main() {
    printf("== 上凸壳.check:与「翻转 y 求下凸壳」「凸包截上半圈」两条参考对拍 + 性质断言 ==\n");

    // ===== 1. 解析用例 =====
    {
        int k;
        vector<p2> sq = {P(0, 0), P(4, 0), P(4, 4), P(0, 4)};
        vector<p2> h = runUpper(sq, 0, k);
        CHECK(k == 2 && h[0] == P(0, 4) && h[1] == P(4, 4), "upper_hull 正方形 -> 上边界两个角(左到右)");
        CHECK(upper_hull(0, sq.data(), sq.data()) == 0, "upper_hull:n=0 -> 0");
        vector<p2> one = {P(3, 7)};
        p2 b1[4];
        CHECK(upper_hull(1, one.data(), b1) == 1 && b1[0] == P(3, 7), "upper_hull:n=1 -> 1 个点");
        vector<p2> two = {P(5, 5), P(1, 1)};
        h = runUpper(two, 0, k);
        CHECK(k == 2 && h[0] == P(1, 1) && h[1] == P(5, 5), "upper_hull:n=2 -> 按 x 排序输出两点");
        {
            vector<p2> samex = {P(2, 1), P(2, 9), P(0, 0), P(4, 0), P(2, 5)};  // 同一 x 上三点
            h = runUpper(samex, 0, k);
            CHECK(k == 3 && h[0] == P(0, 0) && h[1] == P(2, 9) && h[2] == P(4, 0), "upper_hull:同一 x 多点只留最高的那个");
        }
        {
            vector<p2> top;  // 矩形顶边塞满共线点
            For(i, 0, 10) top.push_back(P(i, 10));
            top.push_back(P(0, 0)), top.push_back(P(10, 0)), top.push_back(P(5, 3));
            h = runUpper(top, 0, k);
            CHECK(k == 2 && h[0] == P(0, 10) && h[1] == P(10, 10), "upper_hull nos=0:顶边上 11 个共线点 -> 只剩两端");
            h = runUpper(top, 1, k);
            CHECK(k == 11, "upper_hull nos=1:顶边上 11 个共线点全部保留(k=11)");
            int bad = 0;
            ForD(i, 1, 11) if(!(h[i] == P(i, 10))) ++bad;
            CHECK(bad == 0, "upper_hull nos=1:保留的正好是顶边上的 11 个点、按 x 递增");
        }
        {
            vector<p2> col;  // 全共线
            For(i, -5, 5) col.push_back(P(i, 2 * i + 1));
            h = runUpper(col, 0, k);
            CHECK(k == 2 && h[0] == P(-5, -9) && h[1] == P(5, 11), "upper_hull 全共线点 -> 只剩两端点");
            vector<p2> allsame(20, P(3, 3));
            h = runUpper(allsame, 0, k);
            CHECK(k == 1 && h[0] == P(3, 3), "upper_hull 全同点 -> 只剩一个点");
            vector<p2> ver = {P(1, 1), P(1, 5), P(1, 9)};  // 全在一条竖直线上
            h = runUpper(ver, 0, k);
            CHECK(k == 1 && h[0] == P(1, 9), "upper_hull 全在一条竖直线上的点 -> 只剩最高的那个");
        }
        {
            vector<p2> cir;  // 圆上点:上凸壳 = 上半圆上的点(严格凸)
            ForD(i, 0, 24) cir.push_back(P(1000 * cos(2 * pi * i / 24), 1000 * sin(2 * pi * i / 24)));
            h = runUpper(cir, 0, k);
            int bad = 0;
            if(k != 13) bad = 1;  // 从 (1000,0) 逆时针到 (-1000,0) 共 13 个点(y >= 0)
            ForD(i, 0, k) if(cmp(h[i].y, 0) < 0 && cmp(h[i].y, -1e-6L) < 0) bad = 1;
            CHECK(bad == 0, "upper_hull 半径 1000 圆上 24 点 -> 上半圆的 13 个点(严格凸,无共线)");
        }
    }

    // ===== 2/3/4. 随机对拍 =====
    {
        int bad = 0;
        For(t, 1, 4000) {
            int n = (int) rnd(1, 40), R = (int) rnd(1, 12);
            vector<p2> v;
            ForD(i, 0, n) v.push_back(P((db) rnd(-R, R), (db) rnd(-R, R)));
            if(!oneCase(v, bad < 3) && ++bad >= 3) break;
        }
        CHECK(bad == 0, "upper_hull 与两条独立参考逐位一致 + 全部链性质(4000 组:n=1..40,坐标 ±12,大量重复与共线)");
    }
    {
        int bad = 0;
        For(t, 1, 3000) {
            int n = (int) rnd(1, 60);
            vector<p2> v;
            ForD(i, 0, n) v.push_back(P((db) rnd(-1000000, 1000000), (db) rnd(-1000000, 1000000)));
            if(!oneCase(v, bad < 3) && ++bad >= 3) break;
        }
        CHECK(bad == 0, "upper_hull 与参考一致(3000 组:n=1..60,坐标 ±1e6 一般位置)");
    }
    {
        int bad = 0;  // 上边界共线点很多的构造用例
        For(t, 1, 3000) {
            int m = (int) rnd(1, 20);
            vector<p2> v;
            ForD(i, 0, m + 1) v.push_back(P((db) rnd(-100, 100), 50));      // 一堆共线「顶部」点(同一 y)
            ForD(i, 0, m) v.push_back(P((db) rnd(-100, 100), rnd(-100, 20)));  // 下方随机点
            ForD(i, 0, 3) v.push_back(P((db) rnd(-100, 100), 50));           // 再补几个顶点
            int k0 = 0, k1 = 0;
            vector<p2> h0 = runUpper(v, 0, k0), h1 = runUpper(v, 1, k1);
            db mx = 100;
            if(!chainProps(v, h0, 0, mx) || !chainProps(v, h1, 1, mx)) {
                if(bad < 3) prv("输入", v), prv("nos=0", h0), prv("nos=1", h1);
                ++bad;
                continue;
            }
            if(!sameChain(h0, refUpper(v))) {  // nos=0 与参考一致(nos=1 的共线保留语义单独判)
                if(bad < 3) printf("      [反例] nos=0 与参考不同\n"), prv("输入", v), prv("nos=0", h0), prv("参考", refUpper(v));
                ++bad;
                continue;
            }
            // nos=1:所有「y 最高且是该 x 上最高点」的点都要在链上(重复点/同 x 只算一个)
            db mxy = -1e30L;
            ForD(i, 0, (int) v.size()) mxy = max(mxy, v[i].y);
            vector<db> tx;
            ForD(i, 0, (int) v.size()) if(!cmp(v[i].y, mxy)) tx.push_back(v[i].x);
            sort(tx.begin(), tx.end());
            tx.erase(unique(tx.begin(), tx.end()), tx.end());
            int cntTop = tx.size();
            if(k1 < cntTop) {
                if(bad < 3) printf("      [反例] nos=1 的链(%d 点)少于最高一层的不同 x 个数(%d)\n", k1, cntTop);
                ++bad;
            }
        }
        CHECK(bad == 0, "upper_hull:上边界共线点很多的构造用例(3000 组,每组 1..20 个共线顶点)—— 链性质 + 参考一致 + nos=1 保留全部顶层点");
    }

    // ===== 5. ±1e9 =====
    {
        int bad = 0;
        For(t, 1, 1000) {
            int n = (int) rnd(1, 50);
            vector<p2> v;
            ForD(i, 0, n) v.push_back(P((db) rnd(-1000000000, 1000000000), (db) rnd(-1000000000, 1000000000)));
            if(!oneCase(v, bad < 3) && ++bad >= 3) break;
        }
        CHECK(bad == 0, "upper_hull:±1e9 大整数坐标与参考逐位一致(1000 组;整数叉积精确)");
    }
    PASSED("上凸壳");
}
