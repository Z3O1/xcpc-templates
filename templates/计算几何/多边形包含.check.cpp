// 多边形包含 自测:wn_contain(缠绕数非零规则)与「精确整数缠绕数」「有向角积分」「射线法奇偶」四方对照
//
// 参考实现(全部自己写,不用模板):
//   exactWN:整数坐标下用 __int128 精确算缠绕数与边界(叉积 == 0 且落在包围盒内 -> 边界),
//     给出 0/1/2 的**真值**(对自交多边形同样成立,缠绕数定义与实现无关)。
//   angWN:把每条边相对 p 的有向角(atan2 差,归一到 (-π,π])累加 = 2πk,k 就是缠绕数 ——
//     与「跨 y 符号 + 叉积符号」的实现完全不同的算法,用来做第二意见。
//   geo.cpp 的 contain(射线法,奇偶规则):简单多边形上必须与本模板一致;自交多边形上
//     两者**故意不同**(非零规则 vs 奇偶规则),这一段只打印观察不改判定。
//
// 用例规模:
//   1 解析用例:正方形网格扫描、凹多边形(梳子形)、蝴蝶结的两个翅膀、顶点/边上点、
//     相邻重复点(退化边)、全共线退化多边形、n<3
//   2 3000 组随机简单多边形(星形,逆/顺时针各半)+ 每组 80 个随机整数查询点,与 exactWN/angWN/contain 四方对照
//   3 3000 组随机**可能自交**的多边形(随机点随机顺序)+ 每组 80 个随机整数查询点,与 exactWN/angWN 对照
//   4 3000 组 ±1e9 大整数坐标(简单的与自交的都有),与 exactWN 对照(整数叉积精确,不受 eps 影响)
//   5 性质级:平移不变性、逆序(反向遍历)不改答案、顶点顺序循环移位不改答案
#include "../_check_base.hpp"
#include "geo.cpp"
#include "多边形包含.cpp"

static p2 P(db x, db y) { return {x, y}; }
static db Abs(db x) { return x < 0 ? -x : x; }
static void pr(p2 a) { printf("(%g,%g)", (double) a.x, (double) a.y); }
static void prv(const char *s, const vector<p2> &v) {
    printf("      %s(%zu):", s, v.size());
    ForD(i, 0, (int) min<size_t>(v.size(), 30)) printf(" "), pr(v[i]);
    printf("%s\n", v.size() > 30 ? " ..." : "");
}

// ---------- 独立参考实现 ----------
static int exactWN(const vector<p2> &a, ll px, ll py) {  // 整数坐标下的精确缠绕数 + 边界:0/1/2
    int n = a.size(), wn = 0;
    ForD(i, 0, n) {
        ll ux = (ll) a[i].x, uy = (ll) a[i].y, vx = (ll) a[(i + 1) % n].x, vy = (ll) a[(i + 1) % n].y;
        if(ux == vx && uy == vy) continue;
        i128 cr = (i128) (vx - ux) * (py - uy) - (i128) (vy - uy) * (px - ux);
        if(cr == 0 && min(ux, vx) <= px && px <= max(ux, vx) && min(uy, vy) <= py && py <= max(uy, vy)) return 1;
        if(uy <= py) {
            if(vy > py && cr > 0) ++wn;
        } else if(vy <= py && cr < 0)
            --wn;
    }
    return wn ? 2 : 0;
}
static bool angWNOk(const vector<p2> &a, p2 p, int &sign_out) {  // 有向角和法;false 表示退化(p 在边上/顶点上)
    int n = a.size();
    db s = 0;
    ForD(i, 0, n) {
        p2 u = a[i] - p, v = a[(i + 1) % n] - p;
        if(dis2(u) < 1e-18L || dis2(v) < 1e-18L) return false;
        db d = remainderl(v.alpha() - u.alpha(), 2 * pi);
        if(Abs(Abs(d) - pi) < 1e-6L) return false;  // 边穿过 p 所在直线:可能贴着边界,排除
        s += d;
    }
    db k = s / (2 * pi), kr = roundl(k);  // k 是缠绕数;s 会有 1e-20 量级的舍入残差,必须先 round 再定符号
    if(Abs(k - kr) > 1e-6L) return false;
    sign_out = kr > 0 ? 1 : kr < 0 ? -1 : 0;
    return true;
}
static bool onAnyEdge(const vector<p2> &a, p2 p) {  // 与模板同约定:跳过退化的零长边(否则 geo 的 ons 对任意点都为真)
    int n = a.size();
    ForD(i, 0, n) if(!(a[i] == a[(i + 1) % n]) && ons({a[i], a[(i + 1) % n]}, p)) return true;
    return false;
}

int main() {
    printf("== 多边形包含.check:缠绕数 vs 精确整数缠绕数/有向角积分/射线法 ==\n");

    // ===== 1. 解析用例 =====
    {
        vector<p2> sq = {P(0, 0), P(4, 0), P(4, 4), P(0, 4)};
        int bad = 0;
        For(x, -2, 6) For(y, -2, 6) {
            p2 p = P((db) x / 2, (db) y / 2);
            int got = wn_contain(4, sq.data(), p);
            bool in = x >= 0 && x <= 8 && y >= 0 && y <= 8;
            bool bd = in && (x == 0 || x == 8 || y == 0 || y == 8);
            int want = bd ? 1 : (in ? 2 : 0);
            if(got != want) {
                if(!bad) printf("      [反例] 正方形网格 "), pr(p), printf(":模板 %d,期望 %d\n", got, want);
                ++bad;
            }
        }
        CHECK(bad == 0, "wn_contain 正方形 9×9 网格:0(外)/1(边与顶点)/2(内)全部正确");
        vector<p2> cw = {P(0, 0), P(0, 4), P(4, 4), P(4, 0)};  // 顺时针
        CHECK(wn_contain(4, cw.data(), P(2, 2)) == 2 && wn_contain(4, cw.data(), P(2, 0)) == 1 && wn_contain(4, cw.data(), P(5, 2)) == 0,
              "wn_contain 与顶点顺序无关(顺时针多边形同样正确)");
        // 凹多边形(梳子形):逐点与精确整数缠绕数对照
        vector<p2> comb = {P(0, 0), P(6, 0), P(6, 1), P(5, 1), P(5, 5), P(4, 5), P(4, 1), P(3, 1), P(3, 5), P(2, 5), P(2, 1), P(1, 1), P(1, 5), P(0, 5)};
        bad = 0;
        For(x, -1, 7) For(y, -1, 6) {
            p2 p = P(x, y);
            int got = wn_contain(comb.size(), comb.data(), p), want = exactWN(comb, x, y);
            if(got != want) {
                if(!bad) printf("      [反例] 梳子形 "), pr(p), printf(":模板 %d,精确 %d\n", got, want);
                ++bad;
            }
        }
        CHECK(bad == 0, "wn_contain 凹多边形(梳子形)8×7 网格逐点与精确缠绕数一致(含凹口/边界/角点)");
        // 蝴蝶结(自交):两个翅膀的三角形内部在 (1,2)/(3,2),交叉点 (2,2) 在两条边上
        vector<p2> bow = {P(0, 0), P(4, 4), P(4, 0), P(0, 4)};
        CHECK(wn_contain(4, bow.data(), P(1, 2)) == 2 && wn_contain(4, bow.data(), P(3, 2)) == 2,
              "wn_contain 自交蝴蝶结:两个翅膀内部都判 2(缠绕数非零)");
        CHECK(wn_contain(4, bow.data(), P(2, 2)) == 1 && wn_contain(4, bow.data(), P(1, 1)) == 1 && wn_contain(4, bow.data(), P(1, 3)) == 1,
              "wn_contain 自交蝴蝶结:交叉点 (2,2) 与两条对角线上的点 -> 边界 1");
        {
            // 「同一个三角形绕两圈」:缠绕数 2(非零规则判内部),奇偶规则判外部 —— 两种规则的分水岭
            vector<p2> dbl = {P(0, 0), P(4, 0), P(0, 4), P(0, 0), P(4, 0), P(0, 4)};
            CHECK(wn_contain(6, dbl.data(), P(1, 1)) == 2 && exactWN(dbl, 1, 1) == 2,
                  "wn_contain 绕两圈的三角形:缠绕数 2 -> 判内部(与精确缠绕数一致)");
            CHECK(contain(6, dbl.data(), P(1, 1)) == 0, "同一个点:geo.cpp 的 contain(奇偶规则)判外部 —— 两种规则确实不同");
        }
        {
            // 五角星(每两个顶点连一条):中心小五边形缠绕数 ±2 -> 非零规则判内部,奇偶规则判外部
            vector<p2> st;
            ForD(i, 0, 5) st.push_back(P(100 * cos(pi / 2 + 2 * pi * i / 5), 100 * sin(pi / 2 + 2 * pi * i / 5)));
            vector<p2> pent = {st[0], st[2], st[4], st[1], st[3]};
            int sg = 0;
            CHECK(wn_contain(5, pent.data(), P(0, 0)) == 2 && contain(5, pent.data(), P(0, 0)) == 0 && angWNOk(pent, P(0, 0), sg) && sg != 0,
                  "wn_contain 五角星中心:缠绕数 2 判内部,奇偶规则判外部,有向角和法也确认非零");
        }
        // 退化输入
        vector<p2> dupv = {P(0, 0), P(0, 0), P(4, 0), P(4, 4), P(0, 4)};
        CHECK(wn_contain(5, dupv.data(), P(2, 2)) == 2 && wn_contain(5, dupv.data(), P(9, 9)) == 0 && wn_contain(5, dupv.data(), P(2, 0)) == 1,
              "wn_contain 相邻重复点(退化边)被跳过:内部/外部/边界判定照常");
        vector<p2> deg = {P(0, 0), P(1, 1), P(2, 2)};
        CHECK(wn_contain(3, deg.data(), P(1, 1)) == 1 && wn_contain(3, deg.data(), P(1, 0)) == 0, "wn_contain 全共线退化多边形:线上点算子边界、线外算外部");
        CHECK(wn_contain(2, sq.data(), P(0, 0)) == 0 && wn_contain(0, sq.data(), P(0, 0)) == 0, "wn_contain:n<3 -> 0");
    }

    // ===== 2. 随机简单多边形:四方对照 =====
    {
        int bad = 0, cnt = 0;
        For(t, 1, 3000) {
            int n = (int) rnd(3, 12);
            db R = (db) rnd(3, 40);
            vector<p2> a;
            ForD(i, 0, n) {  // 以原点为心按角度排序 -> 简单多边形(星形),再取整
                db ang = 2 * pi * i / n + (db) rnd(-100, 100) / 400 / n;
                db r = R * (db) rnd(40, 100) / 100;
                a.push_back(P(llroundl(r * cos(ang)), llroundl(r * sin(ang))));
            }
            if(t % 2) reverse(a.begin(), a.end());  // 一半顺时针
            bool simple = true;                     // 取整后可能自交/重边,先粗筛:与 geo 的 contain 只对简单多边形比
            ForD(i, 0, n) ForD(j, i + 1, n) if(a[i] == a[j]) simple = false;
            ForD(it, 0, 80) {
                ll qx = rnd(-50, 50), qy = rnd(-50, 50);
                p2 q = P((db) qx, (db) qy);
                int got = wn_contain(n, a.data(), q), exact = exactWN(a, qx, qy);
                ++cnt;
                if(got != exact) {
                    if(bad < 3) printf("      [反例] "), pr(q), printf(":模板 %d,精确 %d\n", got, exact), prv("多边形", a);
                    ++bad;
                    continue;
                }
                int sg = 0;
                if(exact != 1 && angWNOk(a, q, sg)) {  // 第二意见:有向角和法(只对非边界点)
                    int want = sg ? 2 : 0;
                    if(got != want) {
                        if(bad < 3) printf("      [反例] "), pr(q), printf(":模板 %d,有向角和 %d\n", got, want), prv("多边形", a);
                        ++bad;
                    }
                }
                if(simple && exact != 1 && contain(n, a.data(), q) != got) {
                    if(bad < 3) printf("      [反例] "), pr(q), printf(":模板 %d,contain(奇偶) %d(简单多边形上必须一致)\n", got, contain(n, a.data(), q));
                    ++bad;
                }
            }
        }
        CHECK(bad == 0, "wn_contain 与精确缠绕数/有向角和/contain 一致(3000 组随机简单多边形 × 80 点 = 24 万次)");
        printf("  [info] 本段查询次数:%d\n", cnt);
    }

    // ===== 3. 随机可能自交的多边形:与精确缠绕数 + 有向角和对照 =====
    {
        int bad = 0, cnt = 0, inside = 0;
        For(t, 1, 3000) {
            int n = (int) rnd(3, 10);
            vector<p2> a;
            ForD(i, 0, n) a.push_back(P((db) rnd(-15, 15), (db) rnd(-15, 15)));
            ForD(it, 0, 80) {
                ll qx = rnd(-20, 20), qy = rnd(-20, 20);
                p2 q = P((db) qx, (db) qy);
                int got = wn_contain(n, a.data(), q), exact = exactWN(a, qx, qy);
                ++cnt, inside += exact == 2;
                if(got != exact) {
                    if(bad < 3) printf("      [反例] "), pr(q), printf(":模板 %d,精确 %d\n", got, exact), prv("多边形", a);
                    ++bad;
                    continue;
                }
                int sg = 0;
                if(exact != 1 && angWNOk(a, q, sg) && got != (sg ? 2 : 0)) {
                    if(bad < 3) printf("      [反例] "), pr(q), printf(":模板 %d,有向角和 %d\n", got, sg ? 2 : 0), prv("多边形", a);
                    ++bad;
                }
            }
        }
        CHECK(bad == 0, "wn_contain 与精确缠绕数一致(3000 组随机可能自交多边形 × 80 点 = 24 万次)");
        printf("  [info] 本段查询次数:%d(其中非零规则下算内部的 %d 个)\n", cnt, inside);
    }

    // ===== 4. ±1e9 大整数坐标 =====
    {
        int bad = 0;
        For(t, 1, 3000) {
            int n = (int) rnd(3, 10);
            vector<p2> a;
            ForD(i, 0, n) a.push_back(P((db) rnd(-1000000000, 1000000000), (db) rnd(-1000000000, 1000000000)));
            ForD(it, 0, 20) {
                ll qx = rnd(-1000000000, 1000000000), qy = rnd(-1000000000, 1000000000);
                int got = wn_contain(n, a.data(), P((db) qx, (db) qy)), exact = exactWN(a, qx, qy);
                if(got != exact) {
                    if(bad < 3) printf("      [反例] 大坐标 "), pr(P((db) qx, (db) qy)), printf(":模板 %d,精确 %d\n", got, exact), prv("多边形", a);
                    ++bad;
                }
            }
        }
        CHECK(bad == 0, "wn_contain:±1e9 大整数坐标与精确缠绕数一致(3000 组 × 20 点 = 6 万次;整数叉积精确)");
    }

    // ===== 5. 平移/反向/循环移位不变性 =====
    {
        int bad = 0;
        For(t, 1, 1500) {
            int n = (int) rnd(3, 12);
            vector<p2> a;
            ForD(i, 0, n) a.push_back(P((db) rnd(-30, 30), (db) rnd(-30, 30)));
            p2 sh = P((db) rnd(-1000, 1000), (db) rnd(-1000, 1000));
            vector<p2> mv(n), rv(n), cy(n);
            ForD(i, 0, n) mv[i] = a[i] + sh;
            ForD(i, 0, n) rv[i] = a[n - 1 - i];
            ForD(i, 0, n) cy[i] = a[(i + 3) % n];
            ForD(it, 0, 40) {
                p2 q = P((db) rnd(-40, 40), (db) rnd(-40, 40));
                int base = wn_contain(n, a.data(), q);
                if(wn_contain(n, mv.data(), q + sh) != base) ++bad;
                if(wn_contain(n, rv.data(), q) != base) ++bad;
                if(wn_contain(n, cy.data(), q) != base) ++bad;
                bool bd = onAnyEdge(a, q);
                if(bd != (base == 1)) ++bad;  // 边界判定与「确实在某条边上」等价
            }
        }
        CHECK(bad == 0, "wn_contain 平移/反向遍历/循环移位不改结论,且「返回 1」等价于「确实落在某条边上」(1500 组 × 40 点)");
    }
    PASSED("多边形包含");
}
