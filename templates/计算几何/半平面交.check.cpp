// 半平面交 自测:hpi(有向直线左侧求交)与「逐半平面暴力裁剪多边形」独立实现对拍 + 性质断言
//
// 对拍方式(全部为文件内自写的独立参考实现,不用模板):
//   refClip/refHPI:把一个大正方形(-B..B,B=1e4)按每条直线「保留 cross(l,·) >= 0 一侧」逐条裁剪,
//     交点用参数插值 p = u + (v-u)*cu/(cu-cv) 现算(与模板的 isll 加权平均是两套公式)。
//     hpi 的输入统一**补上 B 正方形的 4 条有向边**,于是「无界」也被裁成有界,两侧算的是同一个集合,
//     可以直接逐顶点 + 面积对比。
//   refHull:Andrew 单调链(独立实现),用凸包的每条逆时针边当半平面,结果必须回到凸包本身。
//
// 用例规模:
//   1 基本性质(排序/去平行/任意输入都不崩)与解析用例:正方形、三角形、空交集、退化点/线段、
//     重复直线、同向平行堆叠、反向重合直线……
//   2 2000 组「极角均匀间隔带抖动」的 3~10 条直线 + B 方框,与暴力裁剪逐顶点(1e-9)+ 面积(1e-9 相对)对拍
//   3 2000 组「完全随机方向」直线(跳过相邻极角差 < 0.05 的病态组)对拍
//   4 1500 组随机凸包:用凸包的逆时针边求交,面积/顶点必须与凸包逐位一致
//   5 2000 组「大量重复/平行/冗余直线」的随机对拍(检查去重逻辑)
//   6 所有随机组顺带做性质级断言:顶点数 <= 直线数、逆时针(带号面积为正)、每个顶点都满足每条
//     直线(容差内)、面积不超过方框面积、输出无重合点
//   7 ±1e9 大坐标只做性质级断言!
#include "../_check_base.hpp"
#include "geo.cpp"
#include "半平面交.cpp"

static p2 P(db x, db y) { return {x, y}; }
static db Abs(db x) { return x < 0 ? -x : x; }
static bool eqd(db a, db b, db tol = 1e-9) { return Abs(a - b) < tol; }
static bool eqp(p2 a, p2 b, db tol = 1e-9) { return eqd(a.x, b.x, tol) && eqd(a.y, b.y, tol); }
static void pr(p2 a) { printf("(%g,%g)", (double) a.x, (double) a.y); }
static void prp(const char *s, const vector<p2> &v) {
    printf("      %s(%zu):", s, v.size());
    for(p2 x : v) printf(" "), pr(x);
    printf("\n");
}
static void prl(const char *s, const vector<seg> &v) {
    printf("      %s(%zu):", s, v.size());
    for(seg x : v) printf(" "), pr(x.x), printf("->"), pr(x.y);
    printf("\n");
}
static db polyArea(const vector<p2> &v) {
    db s = 0;
    ForD(i, 0, (int) v.size()) s += v[i].det(v[(i + 1) % v.size()]);
    return s / 2;
}
static p2 gg = P(0, 0);  // 打印反例时固定用的输入缓存(见 runCase)

// ---------- 独立参考实现 ----------
static vector<p2> refClip(vector<p2> poly, seg l) {  // 保留 cross(l, .) >= 0 一侧
    vector<p2> r;
    int n = poly.size();
    if(!n) return r;
    ForD(i, 0, n) {
        p2 u = poly[i], v = poly[(i + 1) % n];
        db cu = cross(l, u), cv = cross(l, v);
        if(cu >= 0) r.push_back(u);
        if((cu > 0 && cv < 0) || (cu < 0 && cv > 0)) r.push_back(u + (v - u) * (cu / (cu - cv)));
    }
    return r;
}
static vector<seg> boxLines(db B) {  // 逆时针正方形 4 条有向边:左侧即内部
    vector<p2> c = {P(-B, -B), P(B, -B), P(B, B), P(-B, B)};
    vector<seg> r;
    ForD(i, 0, 4) r.push_back({c[i], c[(i + 1) % 4]});
    return r;
}
static vector<p2> refHPI(const vector<seg> &ls, db B) {  // 只裁剪 B 方框(调用方自己把方框边并入 ls)
    vector<p2> cur = {P(-B, -B), P(B, -B), P(B, B), P(-B, B)};
    ForD(i, 0, (int) ls.size()) cur = refClip(cur, ls[i]);
    if(cur.size() < 3) return {};
    return cur;
}
static int refCmpP(p2 a, p2 b) {
    if(a.x != b.x) return a.x < b.x ? -1 : 1;
    if(a.y != b.y) return a.y < b.y ? -1 : 1;
    return 0;
}
static vector<p2> refHull(vector<p2> v) {  // Andrew 单调链:去重 + 去共线,逆时针,首点为最小点
    sort(v.begin(), v.end(), [](p2 a, p2 b) { return refCmpP(a, b) < 0; });
    v.erase(unique(v.begin(), v.end(), [](p2 a, p2 b) { return refCmpP(a, b) == 0; }), v.end());
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

// ---------- 规范化:去掉重合点与共线中间点,再按 (x,y) 排序,便于逐元素比较 ----------
static db nearTol(const vector<p2> &v) {  // 与坐标规模相关的容差(仅用于规范化/比较)
    db s = 1;
    for(p2 x : v) s = max(s, max(Abs(x.x), Abs(x.y)));
    return 1e-9 * s;
}
static vector<p2> canon(vector<p2> v) {
    db tol = nearTol(v);
    vector<p2> r;
    ForD(i, 0, (int) v.size()) if(r.empty() || !eqp(r.back(), v[i], tol)) r.push_back(v[i]);
    while(r.size() > 1 && eqp(r.front(), r.back(), tol)) r.pop_back();
    bool ch = 1;
    while(ch && r.size() >= 3) {  // 去掉共线中间点
        ch = 0;
        db lim = 1e-7 * nearTol(r) * nearTol(r);
        vector<p2> t;
        ForD(i, 0, (int) r.size()) {
            p2 a = r[(i + r.size() - 1) % r.size()], b = r[i], c = r[(i + 1) % r.size()];
            if(Abs((b - a).det(c - a)) < lim) {
                ch = 1;
                continue;
            }
            t.push_back(b);
        }
        if(t.size() < 3) {
            r = t;
            break;
        }
        r = t;
    }
    sort(r.begin(), r.end(), [](p2 a, p2 b) { return refCmpP(a, b) < 0; });
    return r;
}
// 把模板输出与参考输出比一遍;返回 false 表示不一致(并打印现场)
static int cmpOut(const vector<seg> &inp, const vector<p2> &got, const vector<p2> &want, db boxArea) {
    vector<p2> A = canon(got), B = canon(want);
    db tol = 1e-9;
    if(A.size() != B.size()) {
        printf("      顶点数不一致(规范去共线后):模板 %zu,参考 %zu\n", A.size(), B.size());
        prp("模板输出", got), prp("参考输出", want), prl("输入直线", inp);
        return 0;
    }
    ForD(i, 0, (int) A.size()) if(!eqp(A[i], B[i], tol)) {
        printf("      第 %d 个顶点不一致(排序后):模板 ", i), pr(A[i]), printf(",参考 "), pr(B[i]), printf("\n");
        prp("模板输出", got), prp("参考输出", want), prl("输入直线", inp);
        return 0;
    }
    db a0 = polyArea(got), a1 = polyArea(want);
    if(!eqd(a0, a1, 1e-9 * (1 + Abs(a1)))) {
        printf("      面积不一致:模板 %.17Lg,参考 %.17Lg\n", a0, a1);
        prp("模板输出", got), prp("参考输出", want), prl("输入直线", inp);
        return 0;
    }
    if(a0 < -1e-9 * (1 + boxArea)) {
        printf("      输出不是逆时针(带号面积 %.17Lg < 0)\n", a0);
        prp("模板输出", got), prl("输入直线", inp);
        return 0;
    }
    if(a0 > boxArea + 1e-6 * (1 + boxArea)) {
        printf("      输出面积 %.17Lg 超过方框面积 %.17Lg\n", a0, boxArea);
        prp("模板输出", got), prl("输入直线", inp);
        return 0;
    }
    ForD(i, 0, (int) got.size()) ForD(j, 0, (int) inp.size()) {  // 每个顶点都要满足每条直线
        db c = cross(inp[j], got[i]);
        if(c < -1e-6 * (1 + Abs(c)) * (1 + dis(inp[j].dir())) * (1 + dis(got[i] - inp[j].x))) {
            printf("      顶点 "), pr(got[i]), printf(" 违反直线 "), pr(inp[j].x), printf("->"), pr(inp[j].y), printf("  cross = %.17Lg\n", c);
            prp("模板输出", got), prl("输入直线", inp);
            return 0;
        }
    }
    return 1;
}

// ---------- 用例 1:解析用例 ----------
static void analytic() {
    {
        vector<seg> sq = {{P(0, 0), P(4, 0)}, {P(4, 0), P(4, 4)}, {P(4, 4), P(0, 4)}, {P(0, 4), P(0, 0)}};
        vector<p2> r = hpi(sq);
        CHECK(r.size() == 4 && eqd(polyArea(r), 16) && eqd(r[0].x, 0) && eqd(r[0].y, 0) && eqd(r[2].x, 4) && eqd(r[2].y, 4),
              "hpi:4 条边围成正方形 -> 4 个顶点、逆时针、面积 16、首点 (0,0)");
        vector<seg> tri = {{P(0, 0), P(6, 0)}, {P(6, 0), P(0, 3)}, {P(0, 3), P(0, 0)}};
        r = hpi(tri);
        CHECK(r.size() == 3 && eqd(polyArea(r), 9), "hpi:三角形 -> 3 个顶点、面积 9");
        vector<seg> red = sq;
        ForD(i, 0, 3) red.push_back(sq[1]);  // 同一条边重复 4 次
        red.push_back({P(4, 0), P(4, 8)});   // 同向、更弱的一条(会切掉正方形?不会:它更靠右)
        r = hpi(red);
        CHECK(r.size() == 4 && eqd(polyArea(r), 16), "hpi:重复直线 + 同向更弱直线 -> 结果不变(去平行只留最强)");
        vector<seg> str = sq;
        str.push_back({P(1, 8), P(1, 0)});  // 与左边 {P(0,4),P(0,0)} 同向的加强版:x >= 1
        r = hpi(str);
        CHECK(r.size() == 4 && eqd(polyArea(r), 12), "hpi:同向更强直线把正方形切掉一条 -> 面积 12");
        vector<seg> emp = sq;
        emp.push_back({P(0, 0), P(0, 1)});  // 左侧是 x <= 0,与正方形 x >= 0 只交于一条边
        emp.push_back({P(1, 0), P(1, 1)});  // 左侧是 x >= 1
        r = hpi(emp);
        CHECK(r.empty(), "hpi:互相矛盾的两条直线 -> 空");
        vector<seg> pt = {seg{P(0, -1), P(0, 1)}, seg{P(0, 1), P(0, -1)}, seg{P(-1, 0), P(1, 0)}, seg{P(1, 0), P(-1, 0)}};
        r = hpi(pt);
        CHECK(r.empty(), "hpi:四条直线只交于一点 -> 退化返回空");
        vector<seg> sgm = {seg{P(0, -1), P(0, 1)}, {P(0, 1), P(0, -1)}, {P(-2, 0), P(2, 0)}};
        r = hpi(sgm);
        CHECK(r.empty(), "hpi:只交于一条线段(x = 0, -2 <= y <= 2)-> 退化返回空");
        r = hpi({});
        CHECK(r.empty(), "hpi:空输入 -> 空");
        vector<p2> b4 = hpi(boxLines(10));
        CHECK(b4.size() == 4 && eqd(polyArea(b4), 400), "hpi:4 条方框边 -> 方框本身(面积 400)");
        vector<seg> bx = boxLines(1);
        bx.push_back({P(-1, -1), P(1, 1)});  // 左下方 -> 切掉一半
        r = hpi(bx);
        CHECK(r.size() == 3 && eqd(polyArea(r), 2), "hpi:方框被对角线切掉一半 -> 三角形、面积 2");
    }
}

// ---------- 用例 2/3/5:与暴力裁剪随机对拍 ----------
static int runRandom(int groups, int mode) {  // mode 0:极角均匀 1:完全随机 2:重复/平行冗余
    const db B = 1e4;
    vector<seg> box = boxLines(B);
    db boxArea = 4 * B * B;
    int bad = 0;
    For(t, 1, groups) {
        vector<seg> ls;
        int k = (int) rnd(3, 9);
        if(mode == 2) k = (int) rnd(3, 12);
        db base = (db) rnd(0, 100000) / 100000 * 2 * pi;
        ForD(i, 0, k) {
            p2 o = P((db) rnd(-10000, 10000) / 100, (db) rnd(-10000, 10000) / 100);
            db ang;
            if(mode == 0)
                ang = base + 2 * pi * i / k + (db) rnd(-1000, 1000) / 100000;
            else if(mode == 1)
                ang = (db) rnd(-3141592, 3141592) / 1000000;
            else
                ang = base + 2 * pi * (i % 3) / 3 + (db) rnd(-200, 200) / 100000;
            ls.push_back({o, o + P(cos(ang), sin(ang))});
        }
        if(mode == 2) {  // 再塞一批重复的、同向平移的、反向的
            ForD(i, 0, 3) {
                seg s = ls[rnd(0, ls.size() - 1)];
                int r = (int) rnd(0, 3);
                if(r == 0) ls.push_back(s);
                else if(r == 1) ls.push_back({s.x + P((db) rnd(-300, 300) / 100, (db) rnd(-300, 300) / 100), s.y + P((db) rnd(-300, 300) / 100, (db) rnd(-300, 300) / 100)});
                else if(r == 2) ls.push_back({s.y, s.x});  // 反向 -> 交集成一条线,退化
            }
        }
        if(mode == 1) {  // 病态组(相邻极角差过小)直接跳过,顺便断言模板不崩
            vector<db> a;
            for(seg s : ls) a.push_back(s.dir().alpha());
            sort(a.begin(), a.end());
            bool ok = 1;
            ForD(i, 0, (int) a.size()) if(min(Abs(a[(i + 1) % a.size()] - a[i]), Abs(a[(i + 1) % a.size()] - a[i] - 2 * pi)) < 0.05) ok = 0;
            if(!ok) {
                hpi(ls);
                continue;
            }
        }
        vector<seg> inp = ls;
        ForD(i, 0, 4) inp.push_back(box[i]);
        vector<p2> got = hpi(inp), want = refHPI(ls, B);
        if(!cmpOut(inp, got, want, boxArea)) {
            printf("      [反例] mode=%d 第 %d 组(共 %d 条直线,未含方框)\n", mode, t, (int) ls.size());
            if(++bad >= 3) break;
        }
    }
    return bad;
}

// ---------- 用例 4:凸包的每条逆时针边当半平面 -> 必须回到凸包 ----------
static int runHull(int groups) {
    int bad = 0;
    For(t, 1, groups) {
        int n = (int) rnd(3, 40);
        vector<p2> v;
        ForD(i, 0, n) v.push_back(P((db) rnd(-1000, 1000), (db) rnd(-1000, 1000)));
        vector<p2> H = refHull(v);
        if(H.size() < 3) continue;
        vector<seg> ls;
        ForD(i, 0, (int) H.size()) ls.push_back({H[i], H[(i + 1) % H.size()]});
        vector<p2> got = hpi(ls);
        if(!cmpOut(ls, got, H, 4e6)) {
            printf("      [反例] 凸包对拍第 %d 组(凸包 %zu 个顶点)\n", t, H.size());
            prp("输入凸包", H);
            if(++bad >= 3) break;
        }
    }
    return bad;
}

// ---------- 用例 7:±1e9 大坐标只做性质断言 ----------
static int runBig(int groups) {
    int bad = 0;
    For(t, 1, groups) {
        int k = (int) rnd(3, 8);
        vector<seg> ls;
        db base = (db) rnd(0, 100000) / 100000 * 2 * pi;
        ForD(i, 0, k) {
            ll x = rnd(-1000000000, 1000000000), y = rnd(-1000000000, 1000000000);
            db ang = base + 2 * pi * i / k + (db) rnd(-1000, 1000) / 100000;
            p2 o = P((db) x, (db) y);
            ls.push_back({o, o + P(cos(ang) * 1000, sin(ang) * 1000)});
        }
        vector<p2> got = hpi(ls);
        if(got.empty()) continue;
        db a = polyArea(got);
        // 大坐标下 eps = 1e-10 的绝对叉积容差很松,只做「面积非负 + 顶点满足直线(松容差)」两条
        if(a < 0) {
            printf("      大坐标下带号面积为负 %.17Lg\n", a);
            prl("输入", ls), prp("输出", got);
            if(++bad >= 3) break;
            continue;
        }
        ForD(i, 0, (int) got.size()) ForD(j, 0, (int) ls.size()) {
            db c = cross(ls[j], got[i]), sc = (1 + dis(ls[j].dir())) * (1 + dis(got[i] - ls[j].x));
            if(c < -1e-6L * sc) {
                printf("      大坐标下顶点违反直线:cross = %.17Lg,规模 %.3Lg\n", c, sc);
                prl("输入", ls), prp("输出", got);
                if(++bad >= 3) break;
            }
        }
        if(bad >= 3) break;
    }
    return bad;
}

int main() {
    printf("== 半平面交.check:暴力裁剪对拍 + 凸包对拍 + 解析/退化用例 ==\n");
    analytic();

    int b2 = runRandom(2000, 0);
    CHECK(b2 == 0, "hpi 与暴力裁剪一致:(2000 组极角均匀间隔的 3~9 条直线 + B=1e4 方框,逐顶点 1e-9 + 面积 1e-9 相对)");
    int b3 = runRandom(2000, 1);
    CHECK(b3 == 0, "hpi 与暴力裁剪一致:(2000 组随机方向直线,跳过相邻极角差 < 0.05 的病态组)");
    int b5 = runRandom(2000, 2);
    CHECK(b5 == 0, "hpi 与暴力裁剪一致:(2000 组含重复/同向平移/反向重合直线的冗余输入)");

    int b4 = runHull(1500);
    CHECK(b4 == 0, "hpi:把随机凸包的每条逆时针边当半平面 -> 顶点/面积与原凸包逐位一致(1500 组,凸包 3~40 顶点)");

    int b7 = runBig(3000);
    CHECK(b7 == 0, "hpi:±1e9 大坐标下性质成立(面积非负、顶点满足每条直线(松容差),3000 组)");

    {  // 顶点数上界与「无重合点」:任意输入下输出顶点数不超过输入直线数
        int bad = 0;
        For(t, 1, 2000) {
            int k = (int) rnd(3, 12);
            vector<seg> ls;
            ForD(i, 0, k) {
                p2 o = P((db) rnd(-500, 500), (db) rnd(-500, 500));
                db ang = (db) rnd(-3141592, 3141592) / 1000000;
                ls.push_back({o, o + P(cos(ang), sin(ang))});
            }
            vector<p2> got = hpi(ls);
            if((int) got.size() > k) ++bad;
            ForD(i, 0, (int) got.size()) if(eqp(got[i], got[(i + 1) % got.size()], 1e-12)) ++bad;
            if((int) got.size() >= 3 && polyArea(got) <= 0) ++bad;
        }
        CHECK(bad == 0, "hpi 性质:输出顶点数 <= 直线数、无重合点、逆时针(2000 组随机)");
    }
    PASSED("半平面交");
}
