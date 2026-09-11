// geo 自测(强化版):p2 点/向量运算、比较函数、距离、线段位置关系、凸包、圆相关等
//   与内置的独立参考实现 + 暴力对照大规模随机对拍,并做大量「性质级」断言。
//
// ⚠ 本 check 目前是 FAIL 状态(预期),原因是 templates/计算几何/geo.cpp 第 33、34 行:
//     p2 operator-=(p2 &x, p2 y) { return x = x + y; }   // 应为 x = x - y
//     p2 operator/=(p2 &x, p2 y) { return x = x + y; }   // 应为 x = x / y,且第二参数应为 db
//   也就是「减等」与「除等」都写成了加法 —— 本文件末尾「已知 bug 复现」段会先把这两个运算符
//   的实测值原样打印出来,再让断言失败(exit 1)。修好 geo.cpp 第 33/34 行后本文件应变成 PASSED
//   (已用 tmp/ 下的补丁副本实测验证);若仓库主人选择直接删掉 /= 重载,需要同步改本文件最后一段。
//
//   注:第 34 行的签名只有 (p2&, p2),连 (p2&, db) 版都没有 —— 写 `a /= 2.0` 会直接编译失败
//   (error: no match for 'operator/='),所以这里用 SFINAE 探测该重载是否存在,再单独测
//   (p2&, p2) 版的实测行为,以保证 check 自身在修 bug 前后都能编译通过。
//
// 覆盖点(每段给出用例规模;整数坐标处与参考实现逐位对齐,1e9 量级处只做性质级断言):
//   1  eps/pi/sign/cmp
//   2  p2 构造/==/</+-*/ 与手算
//   3  r90/rot(含 3000 组随机:公式、保长、角度加成)
//   4  cross/crossop/inc
//   5  isMid/ons/ons_s
//   6  proj/reflect/nearest(3 万组随机 + 零长线段退化观察)
//   7  chkss/chkss_s/chkll/isll/eqll(3 万组随机解析对照 + 退化观察)
//   8  disss(3 万组随机)
//   9  area/contain:正方形网格、凹多边形独立射线法、顺/逆时针、边界点 0/1/2 语义、退化多边形观察
//  10  convex_hull:Andrew 参考对照(3000 组)、输入点数上界、全部共线/全部重合/含重复点、
//      圆上点、大量共线点、n=0/1/2/3、旋转/平移/缩放不变性、nos=1 语义;convex_diameter(暴搜最远点对)
//  11  convex_cut:参考半平面裁剪对照 + 面积守恒/凸性/旋向保持/顶点只能来自输入或切割线
//  12  circumcircle_diameter/circumcenter:2R 关系、圆心等距、斜边中点、共线返回 -1、最小圆覆盖暴力
//  13  1e9 量级数值抗性:hull/area/contain/nearest/disss/proj/ons/chkss/isll/rot 只做性质断言
//  14  1e9 量级 + 大坐标圆相关
//  15  已知 bug 复现(-= 与 /=),以及 Graham 扫描在 -= 出错时的实际后果
//
// 本次运行结论:+= 正确;其余接口在上面各段 [ok](126 条 [ok] 断言、约 30 万次随机对照,单次跑约 0.25s),
//   只有 -= 与 /= 复现了已知 bug,故恒为 FAIL。已验证本 check 不是「恒失败」:
//   bash tmp/geo_verify_patched.sh —— 把 geo.cpp:33-34 修好后,同一份 check 输出 PASSED geo(exit 0)。
//   另做过变异测试:在补丁副本上分别注入 8 个错误(proj 分母 +1e-6、area 面积 ×1.000001、chkss 端点相接
//   改严格、convex_diameter 结果 ×0.999999、外接圆直径 ×1.000001、convex_cut 去掉 d1 == 0 分支、
//   nearest 距离 ×1.000001、ons 漏掉端点),8 个全部被本 check 抓到并 FAIL。
//
// 另:第 6/7/9/10/11 段末尾有几条 [note](不计入失败),是本轮新发现的「退化输入」行为,已单独记为发现项
//   (详见段内说明),不需要改模板就能复现,不影响本 check 的 FAIL/PASS 判定。要点:
//     · ons/isMid 在零长线段上退化为「任意点都在区间内」;nearest/proj 在零长线段上返回 NaN
//       (proj 里 dis2(dir) == 0 -> 0/0);disss 同样被带成 NaN
//     · chkss 把退化成点的线段当作「与投影重叠的线段相交」,不看 crossop 的真实值
//     · convex_hull(n <= 1) 返回 1 但不写 b;输入点全部重合时返回 k = 2 且两个顶点相同
//     · convex_hull(nos=1) 的「保留共线点」受 crossop 的绝对 eps 限制,大坐标下的严格共线点可能仍被丢掉
//     · convex_cut 的输出点数可达 n+1(正常上界),且切割线端点顺序反过来可能多/少一个等价点
#include "../_check_base.hpp"
#include "geo.cpp"

// ---- 小工具 ----
static p2 P(db x, db y) { return {x, y}; }
static db Abs(db x) { return x < 0 ? -x : x; }
static bool eqd(db a, db b, db tol = 1e-9) { return Abs(a - b) < tol; }
static bool eqp(p2 a, p2 b, db tol = 1e-9) { return eqd(a.x, b.x, tol) && eqd(a.y, b.y, tol); }
static db dd(p2 a, p2 b) { return dis2(a - b); }
static void pr(p2 a) { printf("(%g,%g)", (double) a.x, (double) a.y); }

// 失败时把细节打出来的小工具:bad 计数 + 第一条反例的现场
struct Probe {
    int bad = 0, shown = 0;
    void hit(const char *what, db got, db want) {
        if(!bad) printf("      首个反例:%s:实测 %.17Lg,期望 %.17Lg\n", what, got, want);
        ++bad;
    }
    void hitp(const char *what, p2 got, p2 want) {
        if(!bad) {
            printf("      首个反例:%s:实测 ", what);
            pr(got);
            printf(",期望 ");
            pr(want);
            printf("\n");
        }
        ++bad;
    }
};

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

// ================= 独立参考实现(自写,不依赖模板的凸包/包含/裁剪) =================
// Andrew 单调链:精确整数比较(bits/stdc++ 的 sort,不用模板的 eps 版 operator<)
static int refCmpP(p2 a, p2 b) {
    if(a.x != b.x) return a.x < b.x ? -1 : 1;
    if(a.y != b.y) return a.y < b.y ? -1 : 1;
    return 0;
}
static db refDet(p2 a, p2 b) { return a.x * b.y - a.y * b.x; }
static vector<p2> refHull(vector<p2> v) {  // 去重 + 去共线,逆时针,首位为最小点
    sort(v.begin(), v.end(), [](p2 a, p2 b) { return refCmpP(a, b) < 0; });
    v.erase(unique(v.begin(), v.end(), [](p2 a, p2 b) { return refCmpP(a, b) == 0; }), v.end());
    int n = v.size();
    if(n <= 2) return v;
    vector<p2> h;
    ForD(i, 0, n) {
        while(h.size() > 1 && sign(refDet(h[h.size() - 1] - h[h.size() - 2], v[i] - h[h.size() - 2])) <= 0) h.pop_back();
        h.push_back(v[i]);
    }
    int t = h.size();
    rFor(i, n - 2, 0) {
        while((int) h.size() > t && sign(refDet(h[h.size() - 1] - h[h.size() - 2], v[i] - h[h.size() - 2])) <= 0) h.pop_back();
        h.push_back(v[i]);
    }
    h.pop_back();
    return h;
}
static int refOnSeg(p2 a, p2 b, p2 q) {  // q 在线段 ab 上(整数坐标,det 精确)
    return sign(refDet(b - a, q - a)) == 0 && min(a.x, b.x) <= q.x && q.x <= max(a.x, b.x) && min(a.y, b.y) <= q.y && q.y <= max(a.y, b.y);
}
static int refContain(const vector<p2> &h, p2 q) {  // 单位正方形射线法上取整:0 外 / 1 边界 / 2 内
    int n = h.size();
    ForD(i, 0, n) if(refOnSeg(h[i], h[(i + 1) % n], q)) return 1;
    int c = 0;
    ForD(i, 0, n) {
        p2 u = h[i], v = h[(i + 1) % n];
        if(u.y > v.y) swap(u, v);  // 半开区间 (u.y, v.y]:水平边被排除、顶点只数一次
        if(!(u.y < q.y && q.y <= v.y)) continue;
        if(refDet(v - u, q - u) > 0) ++c;
    }
    return c % 2 ? 2 : 0;
}
static db refArea(const vector<p2> &h) {  // 带号面积(shoelace)
    db s = 0;
    int n = h.size();
    ForD(i, 0, n) s += refDet(h[i], h[(i + 1) % n]);
    return s / 2;
}
static vector<p2> refCut(const vector<p2> &h, seg q) {  // 保留 cross(q, .) >= 0 一侧
    vector<p2> r;
    int n = h.size();
    ForD(i, 0, n) {
        p2 p1 = h[i], p2r = h[(i + 1) % n];
        db c1 = cross(q, p1), c2 = cross(q, p2r);
        if(c1 >= 0) r.push_back(p1);
        if((c1 > 0 && c2 < 0) || (c1 < 0 && c2 > 0)) r.push_back(p1 + (p2r - p1) * (c1 / (c1 - c2)));
    }
    return r;
}
static db refNearest(p2 a, p2 b, p2 q) {  // 点到线段距离平方
    if(a.x == b.x && a.y == b.y) return dd(a, q);
    p2 dir = b - a;
    db t = (dir * (q - a)) / dis2(dir);
    return dd(q, a + dir * max((db) 0, min((db) 1, t)));
}
static bool refChkss(p2 a, p2 b, p2 c, p2 e) {  // 闭线段相交(参数法 + 共线投影精确)
    p2 d = b - a, f = e - c, g = c - a;
    db den = refDet(d, f);
    if(Abs(den) > 1e-12) {
        db t = refDet(g, f) / den, u = refDet(g, d) / den;
        return t >= -1e-12 && t <= 1 + 1e-12 && u >= -1e-12 && u <= 1 + 1e-12;
    }
    if(Abs(refDet(g, d)) > 1e-12) return false;  // 平行不共线
    auto proj = [](p2 p, p2 q, p2 v) { return (p - q) * v; };  // 一维投影坐标
    db p1 = proj(a, c, d), p2v = proj(b, c, d);                // a,b 在方向 d 上的投影
    db q1 = proj(c, c, d), q2 = proj(e, c, d);
    if(p1 > p2v) swap(p1, p2v);
    if(q1 > q2) swap(q1, q2);
    return !(p2v < q1 - 1e-12 || q2 < p1 - 1e-12);
}
static db refDisss(p2 a, p2 b, p2 c, p2 e) {
    if(refChkss(a, b, c, e)) return 0;
    return min({refNearest(a, b, c), refNearest(a, b, e), refNearest(c, e, a), refNearest(c, e, b)});
}
static p2 refProjAnalytic(p2 a, p2 b, p2 q) {  // 参数方程解析解
    db t = ((q - a) * (b - a)) / dis2(b - a);
    return a + (b - a) * t;
}
// 外心解析式(与模板独立的一套写法)
static p2 refCircumcenter(p2 A, p2 B, p2 C) {
    db d = 2 * (A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));
    db ux = ((A.x * A.x + A.y * A.y) * (B.y - C.y) + (B.x * B.x + B.y * B.y) * (C.y - A.y) + (C.x * C.x + C.y * C.y) * (A.y - B.y)) / d;
    db uy = ((A.x * A.x + A.y * A.y) * (C.x - B.x) + (B.x * B.x + B.y * B.y) * (A.x - C.x) + (C.x * C.x + C.y * C.y) * (B.x - A.x)) / d;
    return {ux, uy};
}
// 最小圆覆盖暴力:圆心只能是两点直径中点或三点外心(O(n^3) 枚举,O(n) 验证)—— 用于验证外接圆
static p2 bruteMEC(const vector<p2> &v) {
    int n = v.size();
    if(n == 1) return v[0];
    db best = 1e300L;
    p2 bc = v[0];
    auto tryC = [&](p2 c) {
        db r = 0;
        ForD(i, 0, n) r = max(r, dd(c, v[i]));
        if(r < best) best = r, bc = c;
    };
    ForD(i, 0, n) ForD(j, i + 1, n) tryC((v[i] + v[j]) / 2);
    ForD(i, 0, n) ForD(j, i + 1, n) ForD(k, j + 1, n) {
        if(sign(refDet(v[j] - v[i], v[k] - v[i])) == 0) continue;
        tryC(refCircumcenter(v[i], v[j], v[k]));
    }
    return bc;
}
// 只用模板的 -= 求「相邻点差向量」的单调链下凸壳:用来展示 -= 出 bug 时的实际后果
// (注意:这里刻意不做完整凸包,只取字典序意义下的下凸壳,好让「结果对不对」一目了然)
static int lowerChainWithTemplateMinusEq(int n, p2 *a, p2 *b) {
    sort(a, a + n, [](p2 u, p2 v) { return refCmpP(u, v) < 0; });
    int k = 0;
    ForD(i, 0, n) {
        while(k > 1) {
            p2 rel = b[k - 1];
            rel -= b[k - 2];  // 模板 operator-=:修好后 rel = b[k-1]-b[k-2];出 bug 时 rel = b[k-1]+b[k-2]
            db want = refDet(b[k - 1] - b[k - 2], a[i] - b[k - 2]);
            db got = rel.det(a[i] - b[k - 2]);
            if(got != want) return -1;  // 用模板算出的叉积与真实叉积不同 -> 直接报告错误
            if(want > 0) break;
            --k;
        }
        b[k++] = a[i];
    }
    return k;
}

// ---- 通用断言器 ----
// 随机点集 -> 模板 convex_hull,并按几何性质自检(k <= n、逆时针、凸、首点为最小点、
// 与独立 Andrew 顶序等价、包含所有输入点);返回 false 表示至少一条性质被违反
static bool hullCheckOne(const vector<p2> &v, bool verbose = false) {
    int n = v.size();
    static p2 a[4096], h[4096];
    ForD(i, 0, n) a[i] = v[i];
    int k = convex_hull(n, a, h);
    vector<p2> H(h, h + k), R = refHull(v);
    bool allSame = true;
    ForD(i, 0, n) if(refCmpP(v[i], v[0]) != 0) allSame = false;
    if(n >= 2 && allSame) return k == 2;  // 全同点:模板返回 2 个重复点(第 15 段 [note]),不参与与参考的比对
    if(k > n || k < 0) {
        if(verbose) printf("      k=%d > n=%d\n", k, n);
        return false;
    }
    if(H.empty()) return false;
    p2 mn = v[0];
    ForD(i, 0, n) if(refCmpP(v[i], mn) < 0) mn = v[i];
    if(v.size() == 1) return k == 1;  // n == 1:k 返回 1,但模板不往 b 里写点(见第 10 段 [note]),故只查返回值
    if(!eqp(H[0], mn, 1e-9)) {
        if(verbose) printf("      H[0] 不是最小点\n");
        return false;
    }
    if(n == 2 && !eqp(H[1], v[0] == mn ? v[1] : v[0], 1e-9)) {
        if(verbose) printf("      n=2 两端点不对\n");
        return false;
    }
    if(k >= 3) {  // 逆时针:以 h[0] 为原点看,所有三角形 (h[0],h[i],h[i+1]) 的带号面积都非负
        db sgn2 = 0;
        ForD(i, 0, k) sgn2 += refDet(h[i] - h[0], h[(i + 1) % k] - h[0]);
        db tol = 1e-9L * (1 + Abs(h[0].x) + Abs(h[0].y)) * (1 + Abs(h[1].x - h[0].x) + Abs(h[1].y - h[0].y)) * k;
        if(sgn2 < -tol) {
            if(verbose) printf("      凸包带号面积(以 h[0] 为原点)= %.17Lg < 0,旋向为顺时针\n", sgn2 / 2);
            return false;
        }
    } else {  // k == 2:两个顶点必须不同(全同点输入已在上面单独 return)
        if(refCmpP(h[0], h[1]) == 0) {
            if(verbose) printf("      两个输出顶点相同,但输入点并非全同点\n");
            return false;
        }
    }
    ForD(i, 0, k) {
        p2 u = h[i], w = h[(i + 1) % k], t = h[(i + 2) % k];
        if(refDet(w - u, t - u) < -1e-9L * (1 + Abs(u.x) + Abs(u.y) + Abs(w.x) + Abs(w.y) + Abs(t.x) + Abs(t.y))) {
            if(verbose) printf("      h[%d] 处非凸\n", (i + 1) % k);
            return false;
        }
    }
    ForD(i, 0, n) ForD(j, 0, k) {
        p2 u = h[j], w = h[(j + 1) % k];
        db c = refDet(w - u, v[i] - u);
        if(c < -1e-9L * (1 + Abs(u.x) + Abs(u.y) + Abs(w.x) + Abs(w.y) + Abs(v[i].x) + Abs(v[i].y))) {
            if(verbose) printf("      输入点 "), pr(v[i]), printf(" 在边 h[%d]-h[%d] 外侧\n", j, (j + 1) % k);
            return false;
        }
    }
    // 与独立参考的顶点集合/顺序等价(去共线语义:双方都严格去共线)
    if(k != (int) R.size()) {
        if(verbose) printf("      顶点数 %d != 参考 %d\n", k, (int) R.size());
        return false;
    }
    ForD(i, 0, k) {
        if(refCmpP(h[i], R[i]) != 0) {
            if(verbose) printf("      第 %d 个顶点不同:", i), pr(h[i]), printf(" vs "), pr(R[i]), printf("\n");
            return false;
        }
    }
    return true;
}
// convex_cut 后的多边形面积(用带号面积 -> 绝对值)
static db cutArea(const vector<p2> &v) { return Abs(refArea(v)); }

int main() {
    printf("== geo.check:强化版随机对拍 + 极端/退化用例(独立参考实现见文件头部)==\n");

    // ===== 1. eps / pi / sign / cmp =====
    CHECK(eps == 1e-10, "eps == 1e-10");
    CHECK(eqd(pi, acosl(-1.0L)) && eqd(pi, 3.14159265358979L), "pi == acos(-1)");
    CHECK(sign(0) == 0 && sign(eps / 2) == 0 && sign(-eps / 2) == 0 && sign(eps) == 0, "sign:|x| <= eps 归零");
    CHECK(sign(eps * 10) == 1 && sign(-eps * 10) == -1, "sign:刚超出 eps 取符号");
    CHECK(sign(1e-6L) == 1 && sign(-1e-6L) == -1 && sign(100) == 1 && sign(-100) == -1, "sign:大值取符号");
    CHECK(sign(1e18L) == 1 && sign(-1e18L) == -1 && sign(1e-15L) == 0, "sign:1e18 与 1e-15 边界");
    CHECK(cmp(1, 1) == 0 && cmp(1 + eps / 4, 1) == 0 && cmp(1, 1 + eps / 4) == 0, "cmp:容差内为 0");
    CHECK(cmp(2, 1) == 1 && cmp(1, 2) == -1 && cmp(1, 1 + 1e-6L) == -1, "cmp:容差外取符号");
    {
        int bad = 0;
        For(t, 1, 20000) {
            db x = (db) rnd(-1000000, 1000000) / 1000, y = (db) rnd(-1000000, 1000000) / 1000;
            if(cmp(x, y) != (x < y ? -1 : x > y ? 1 : 0)) ++bad;
            if(sign(x) != (x < 0 ? -1 : x > 0 ? 1 : 0)) ++bad;
        }
        CHECK(bad == 0, "sign/cmp 与精确比较一致(2 万组随机,含整数倍坐标)");
    }

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
    {
        Probe pb;  // 代数恒等式:p2 加法/数乘/点积/叉积的交换/反对称/双线性
        For(t, 1, 20000) {
            p2 x = P((db) rnd(-1000, 1000), (db) rnd(-1000, 1000)), y = P((db) rnd(-1000, 1000), (db) rnd(-1000, 1000)), z = P((db) rnd(-1000, 1000), (db) rnd(-1000, 1000));
            if(!eqp(x + y, y + x)) pb.hit("x+y != y+x", 0, 0);
            if(!eqp((x + y) + z, x + (y + z))) pb.hit("加法不结合", 0, 0);
            if(!eqp(x - y, x + y * -1.0L)) pb.hit("x-y != x+(-y)", 0, 0);
            if(!eqd(x * y, y * x, 1e-9 * (1 + Abs(x * y)))) pb.hit("点积不对称", x * y, y * x);
            if(!eqd(x.det(y), -y.det(x), 1e-9 * (1 + Abs(x.det(y))))) pb.hit("det 不反对称", x.det(y), -y.det(x));
            if(!eqd((x + y).det(z), x.det(z) + y.det(z), 1e-9 * (1 + Abs(x.det(z)) + Abs(y.det(z))))) pb.hit("det 不线性", 0, 0);
            if(!eqp(x * (db) 3 - x * (db) 2, x)) pb.hit("数乘分配", 0, 0);
            if(!eqp((x + y) / 2.0L, x / 2.0L + y / 2.0L)) pb.hit("除以 2 不分配", 0, 0);
            if(!eqd(dis(x), dis(P(0, 0) - x), 1e-9 * (1 + dis(x)))) pb.hit("dis 奇对称", dis(x), dis(P(0, 0) - x));
            if(dis(x) > 1e-6L && !eqd(dis(unit(x)), 1, 1e-9)) pb.hit("unit 后长度不是 1", dis(unit(x)), 1);
            if(!eqd(x * r90(x), 0, 1e-9 * (1 + dis2(x)))) pb.hit("r90 不正交", x * r90(x), 0);
            if(!eqd(dis2(r90(x)), dis2(x), 1e-9 * (1 + dis2(x)))) pb.hit("r90 不保长", dis2(r90(x)), dis2(x));
            if(!eqd(dis2(x), x * x, 1e-9 * (1 + dis2(x)))) pb.hit("dis2 != x*x", dis2(x), x * x);
        }
        CHECK(pb.bad == 0, "p2 代数性质:加法交换/结合、det 反对称与双线性、数乘分配、dis 奇对称、r90 正交保长(2 万组)");
    }

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
    {
        Probe pb;  // 旋转是正交变换:保点积、保 det、保距离、逆旋转可还原
        For(t, 1, 5000) {
            p2 x = P((db) rnd(-1000, 1000), (db) rnd(-1000, 1000)), y = P((db) rnd(-1000, 1000), (db) rnd(-1000, 1000));
            db ang = (db) rnd(-1000000, 1000000) / 100000;
            p2 rx = rot(x, ang), ry = rot(y, ang);
            db sc = 1 + dis2(x) + dis2(y);
            if(!eqd(rx * ry, x * y, 1e-9L * sc)) pb.hit("旋转不保点积", rx * ry, x * y);
            if(!eqd(rx.det(ry), x.det(y), 1e-9L * sc)) pb.hit("旋转不保 det", rx.det(ry), x.det(y));
            if(!eqd(dis(rx - ry), dis(x - y), 1e-9L * (1 + dis(x - y)))) pb.hit("旋转不保距离", dis(rx - ry), dis(x - y));
            if(!eqp(rot(rx, -ang), x, 1e-6L * (1 + dis(x)))) pb.hitp("旋转后再逆旋转不还原", rot(rx, -ang), x);
        }
        CHECK(pb.bad == 0, "rot 是正交变换:保点积/det/距离,逆旋转还原(5000 组)");
    }

    // ===== 4. cross / crossop / inc(左转、右转、共线) =====
    seg L = {P(0, 0), P(2, 0)};
    CHECK(eqd(cross(L, P(1, 1)), 2) && crossop(L, P(1, 1)) == 1, "cross 左转为正");
    CHECK(eqd(cross(L, P(1, -1)), -2) && crossop(L, P(1, -1)) == -1, "cross 右转为负");
    CHECK(crossop(L, P(1, 0)) == 0 && crossop(L, P(5, 0)) == 0 && crossop(L, P(1, eps / 2)) == 0, "cross 共线为 0");
    CHECK(eqd(cross(L, P(1, 1)), L.dir().det(P(1, 1) - L.x)), "cross == dir().det(q - x)");
    CHECK(eqd(cross({P(0, 0), P(3, 4)}, P(0, 0)), 0) && eqd(cross({P(1, 1), P(4, 5)}, P(1, 1)), 0), "cross 起点处必为 0");
    CHECK(inc(L, P(9, 0)) && inc(L, P(0, 0)) && !inc(L, P(0, 1)), "inc 只看直线不看线段范围");
    CHECK(crossop({P(0, 0), P(1, 1)}, P(2, 2)) == 0 && crossop({P(0, 0), P(1, 1)}, P(2, 3)) == 1, "crossop 斜线段左转/共线");
    {
        Probe pb;  // crossop 只是 sign(cross),两者必须处处一致
        For(t, 1, 20000) {
            p2 u = P((db) rnd(-1000, 1000), (db) rnd(-1000, 1000)), v = P((db) rnd(-1000, 1000), (db) rnd(-1000, 1000)), q = P((db) rnd(-1000, 1000), (db) rnd(-1000, 1000));
            if(crossop({u, v}, q) != sign(cross({u, v}, q))) pb.hit("crossop != sign(cross)", crossop({u, v}, q), sign(cross({u, v}, q)));
            if(crossop({u, v}, q) != -crossop({v, u}, q)) pb.hit("crossop 不反对称", crossop({u, v}, q), -crossop({v, u}, q));
            if(cross({u, v}, q) != -cross({v, u}, q)) pb.hit("cross 不反对称", cross({u, v}, q), -cross({v, u}, q));
            if(u.x != v.x && crossop({u, v}, v) != 0) pb.hit("端点 crossop 不为 0", crossop({u, v}, v), 0);
            // 直线上任取参数点的 crossop 必为 0
            db s = (db) rnd(-500, 500) / 100;
            p2 mid = u + (v - u) * s;
            if(crossop({u, v}, mid) != 0 && dis(u - v) > 1) pb.hit("线段上的点 crossop 不为 0", crossop({u, v}, mid), 0);
            // 直线两侧的点 crossop 反号
            p2 nrm = r90(v - u);
            if(dis(u - v) > 1e-6L) {
                p2 o = (u + v) / 2 + unit(nrm) * 10;
                if(crossop({u, v}, o) != 1) pb.hit("法线正侧 crossop 应为 1", crossop({u, v}, o), 1);
                if(crossop({u, v}, (u + v) / 2 - unit(nrm) * 10) != -1) pb.hit("法线负侧 crossop 应为 -1", crossop({u, v}, o), -1);
            }
        }
        CHECK(pb.bad == 0, "cross/crossop 反对称、端点归零、直线两侧反号(2 万组,含线上参数点)");
    }

    // ===== 5. isMid / ons / ons_s =====
    CHECK(isMid(0.0L, 0.5L, 1.0L) && isMid(1.0L, 0.5L, 0.0L), "isMid 区间内(两端顺序无关)");
    CHECK(isMid(0.0L, 0.0L, 1.0L) && isMid(0.0L, 1.0L, 1.0L) && isMid(2.0L, 2.0L, 2.0L), "isMid 端点算在区间内");
    CHECK(!isMid(0.0L, 1.5L, 1.0L) && !isMid(0.0L, -0.5L, 1.0L), "isMid 区间外为假");
    CHECK(isMid(0.0L, 5.0L, 0.0L), "isMid a==b 时任意 m 都为真(退化区间的约定)");
    CHECK(isMid(P(0, 0), P(1, 1), P(2, 2)) && isMid(P(0, 0), P(1, 2), P(2, 2)), "isMid 点版:逐坐标落在区间内(即落在包围盒里)");
    CHECK(!isMid(P(0, 0), P(3, 2), P(2, 2)) && !isMid(P(0, 0), P(-1, 1), P(2, 2)), "isMid 点版:任一坐标出界即为假");
    CHECK(ons(L, P(1, 0)) && ons(L, P(0, 0)) && ons(L, P(2, 0)), "ons 线段内部与端点");
    CHECK(!ons(L, P(3, 0)) && !ons(L, P(-1, 0)) && !ons(L, P(1, eps * 100)), "ons 延长线与线外为假");
    CHECK(ons_s(L, P(1, 0)) && !ons_s(L, P(0, 0)) && !ons_s(L, P(2, 0)) && !ons_s(L, P(1, eps * 100)), "ons_s 严格在线段内部");
    CHECK(ons({P(0, 0), P(0, 3)}, P(0, 1.5)) && ons({P(1, 1), P(3, 3)}, P(2, 2)), "ons 竖直线段/斜线段");
    CHECK(ons({P(-5, -5), P(5, 5)}, P(0, 0)) && !ons({P(-5, -5), P(5, 5)}, P(1, 1 + eps * 1e3)), "ons 单位法向外的点立刻为假");
    {
        Probe pb;  // ons == (crossop == 0 && 在包围盒内);ons_s 再加「严格在端点之间」
        For(t, 1, 30000) {
            p2 u = P((db) rnd(-20, 20), (db) rnd(-20, 20)), v = P((db) rnd(-20, 20), (db) rnd(-20, 20));
            p2 q = P((db) rnd(-25, 25), (db) rnd(-25, 25));
            if(refCmpP(u, v) == 0) continue;  // 零长线段是退化输入,单独观察(见下面 [note])
            int got = ons({u, v}, q);
            int want = (sign(cross({u, v}, q)) == 0) && min(u.x, v.x) <= q.x && q.x <= max(u.x, v.x) && min(u.y, v.y) <= q.y && q.y <= max(u.y, v.y);
            if(got != want) pb.hit("ons 与整数参考不一致", got, want);
            if(ons_s({u, v}, q) != (want && refCmpP(u, q) != 0 && refCmpP(v, q) != 0)) pb.hit("ons_s 与整数参考不一致", ons_s({u, v}, q), 1);
        }
        CHECK(pb.bad == 0, "ons/ons_s 与独立参考一致(3 万组整数坐标,含端点是否计入的严格语义)");
        printf("  [note] ons 依赖的 isMid(p.x, q, p.y) 在零长线段(p.x == p.y)上退化为「任意 q 都在区间内」(不计入失败):"
               "ons({(3,4),(3,4)}, (100,100)) 实测 %d(期望 0);调用方需自行保证线段非退化。\n",
               (int) ons({P(3, 4), P(3, 4)}, P(100, 100)));
    }

    // ===== 6. proj / reflect / nearest(返回平方距离) =====
    seg X = {P(0, 0), P(4, 0)};
    CHECK(eqp(proj(X, P(1, 3)), P(1, 0)) && eqp(proj(X, P(-7, 3)), P(-7, 0)), "proj 垂足(含线段外)");
    CHECK(eqp(reflect(X, P(1, 3)), P(1, -3)) && eqp(reflect(X, P(5, -2)), P(5, 2)), "reflect 镜像");
    CHECK(eqd(nearest(X, P(1, 3)), 9) && eqd(nearest(X, P(-3, 1)), 10) && eqd(nearest(X, P(5, 1)), 2), "nearest 垂足在内/线段外(返回平方距离)");
    CHECK(eqd(nearest(X, P(0, 0)), 0) && eqd(nearest(X, P(4, 0)), 0), "nearest 端点距离为 0");
    CHECK(eqd(nearest({P(0, 0), P(0, 5)}, P(3, 2)), 9), "nearest 竖直线段");
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
    {
        Probe pb;  // proj 的解析式对照 + 垂距最短性 + reflect 是对合且保距
        For(t, 1, 30000) {
            p2 u = P((db) rnd(-1000, 1000), (db) rnd(-1000, 1000)), v = P((db) rnd(-1000, 1000), (db) rnd(-1000, 1000)), q = P((db) rnd(-2000, 2000), (db) rnd(-2000, 2000));
            if(dd(u, v) < 1e-6L) continue;
            seg s = {u, v};
            p2 h = proj(s, q);                                // 垂足(可在线段外)
            p2 ha = refProjAnalytic(u, v, q);                 // 参数方程解析解
            db sc = 1 + dd(h, u);
            if(!eqd(dd(h, ha), 0, 1e-12L * sc)) pb.hit("proj != 解析解", dd(h, ha), 0);
            if(!eqd((q - h) * (v - u), 0, 1e-9L * dis(v - u) * (1 + dis(q - u)))) pb.hit("垂线不垂直于方向", (q - h) * (v - u), 0);
            // 垂足若落在线段内,则它到 q 的距离 <= 端点到 q 的距离(垂距最短)
            if(isMid(u, h, v)) {
                db d1 = dd(h, q), d2 = min(dd(u, q), dd(v, q));
                if(d1 > d2 + 1e-6L * (1 + d2)) pb.hit("垂足不是最近点", d1, d2);
            }
            if(!eqd(dd(reflect(s, q), q), 4 * dd(h, q), 1e-9L * (1 + 4 * dd(h, q)))) pb.hit("reflect 到 q 的距离不是 2*hq", dd(reflect(s, q), q), 4 * dd(h, q));
            if(!eqp(reflect(s, reflect(s, q)), q, 1e-9L * (1 + dis(q - u)))) pb.hitp("reflect 不是对合", reflect(s, reflect(s, q)), q);
            if(!eqd(dis(reflect(s, q) - u), dis(q - u), 1e-9L * (1 + dis(q - u)))) pb.hit("reflect 不保端点距离", dis(reflect(s, q) - u), dis(q - u));
            if(!eqd(nearest(s, q), min(dd(h, q), min(dd(u, q), dd(v, q))) == dd(h, q) && isMid(u, h, v) ? dd(h, q) : min(dd(u, q), dd(v, q)), 1e-9L * (1 + dd(q, u)))) pb.hit("nearest 与垂足/端点公式不一致", nearest(s, q), dd(h, q));
        }
        CHECK(pb.bad == 0, "proj/reflect:解析解一致、垂距最短、reflect 对合保距(3 万组)");
    }
    {  // 退化:零长线段 —— 目前 nearest 返回 NaN(见下面 [note])
        printf("  [note] nearest/proj 在零长线段上返回 NaN(不计入失败):nearest({(3,4),(3,4)}, (0,0)) 实测 %Lg,"
               "此处 proj 的 dis2(dir) == 0 会先算出 0/0。模板调用方请自行保证线段非退化。\n",
               nearest({P(3, 4), P(3, 4)}, P(0, 0)));
        printf("  [note] disss 同样受零长线段影响:disss({(3,4),(3,4)}, {(0,0),(5,12)}) 实测 %Lg(参考值 25 —— 点 (3,4) 到线段的最近距离平方)。\n",
               disss({P(3, 4), P(3, 4)}, {P(0, 0), P(5, 12)}));
    }

    // ===== 7. 线段相交 / 直线交点 =====
    CHECK(chkss({P(0, 0), P(2, 2)}, {P(0, 2), P(2, 0)}) && chkss_s({P(0, 0), P(2, 2)}, {P(0, 2), P(2, 0)}), "chkss/chkss_s 十字相交");
    CHECK(chkss({P(0, 0), P(1, 1)}, {P(1, 1), P(2, 0)}), "chkss 端点相接算相交");
    CHECK(!chkss_s({P(0, 0), P(1, 1)}, {P(1, 1), P(2, 0)}), "chkss_s 端点相接不算规范相交");
    CHECK(chkss({P(0, 0), P(2, 0)}, {P(1, 0), P(3, 0)}) && !chkss({P(0, 0), P(1, 0)}, {P(2, 0), P(3, 0)}), "chkss 共线重叠/共线分离");
    CHECK(!chkss({P(0, 0), P(1, 0)}, {P(2, 0), P(3, 0)}) && !chkss({P(0, 0), P(1, 0)}, {P(0, 2), P(0, 3)}), "chkss 共线但分离 / 平行但分离");
    CHECK(!chkss({P(0, 0), P(1, 0)}, {P(0, 1), P(1, 1)}) && !chkss({P(0, 0), P(1, 1)}, {P(2, 2), P(3, 3)}), "chkss 平行/共线不重叠");
    seg M = {P(0, 1), P(2, 1)}, N = {P(0, 0), P(0, 2)}, L2 = {P(-5, 0), P(5, 0)};
    CHECK(!chkll(L, M) && !chkll(L, L) && chkll(L, N), "chkll 平行(含重合)为 0、相交为 1");
    CHECK(eqll(L, L2) && !eqll(L, M) && !eqll(L, N), "eqll 同一条直线/平行/相交");
    CHECK(eqp(isll(L, N), P(0, 0)), "isll 交点手算 x 轴 ∩ y 轴");
    CHECK(eqp(isll({P(0, 0), P(2, 2)}, {P(0, 2), P(2, 0)}), P(1, 1)), "isll 交点手算 两条对角线");
    CHECK(eqp(isll({P(1, 1), P(3, 1)}, {P(2, -1), P(2, 5)}), P(2, 1)), "isll 交点手算 水平 ∩ 垂直");
    CHECK(eqp(isll({P(0, 0), P(1, 1)}, {P(1, 0), P(0, 1)}), P(0.5L, 0.5L)), "isll 交点落在两条线段的延长线之间(斜率相反)");
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
    {
        Probe pb;  // 相交判定的性质级检查(含端点相接、参数 t/u 落在 [0,1])
        For(t, 1, 20000) {
            p2 u = P((db) rnd(-15, 15), (db) rnd(-15, 15)), v = P((db) rnd(-15, 15), (db) rnd(-15, 15));
            p2 w = P((db) rnd(-15, 15), (db) rnd(-15, 15)), z = P((db) rnd(-15, 15), (db) rnd(-15, 15));
            if(dd(u, v) < 1 || dd(w, z) < 1) continue;
            seg s1 = {u, v}, s2 = {w, z};
            if(chkss_s(s1, s2) && !chkss(s1, s2)) pb.hit("chkss_s 为真但 chkss 为假", 0, 0);
            if(chkss(s2, s1) != chkss(s1, s2)) pb.hit("chkss 不对称", chkss(s1, s2), chkss(s2, s1));
            if(chkss_s(s2, s1) != chkss_s(s1, s2)) pb.hit("chkss_s 不对称", chkss_s(s1, s2), chkss_s(s2, s1));
            // 端点相接必然相交;端点严格在对方线段内部必然相交(ons 语义)
            int touch = ons(s1, w) || ons(s1, z) || ons(s2, u) || ons(s2, v);
            if(touch && !chkss(s1, s2)) pb.hit("有公共点却判不相交", chkss(s1, s2), 1);
            if(ons_s(s1, w) && !chkss(s1, s2)) pb.hit("点严格落在对段内部却判不相交", chkss(s1, s2), 1);
            // chkss_s(规范相交)== 双方严格跨越:某端点严格落在对方线段内部,且对方不含另一个端点
            int d1 = sign(cross(s1, w)), d2 = sign(cross(s1, z));
            int d3 = sign(cross(s2, u)), d4 = sign(cross(s2, v));
            if(chkss_s(s1, s2) != (d1 * d2 < 0 && d3 * d4 < 0)) pb.hit("chkss_s != 端点严格异侧", chkss_s(s1, s2), d1 * d2 < 0 && d3 * d4 < 0);
            if(chkss_s(s1, s2) && !chkss(s1, s2)) pb.hit("规范相交却不被判相交", chkss(s1, s2), 1);
        }
        CHECK(pb.bad == 0, "chkss/chkss_s 对称性与「有公共点即相交」性质(2 万组,含端点相接/内部穿点)");
    }
    {
        Probe pb;  // isll 交出直线上;交点若在线段内则到两段距离为 0
        For(t, 1, 20000) {
            p2 u = P((db) rnd(-20, 20), (db) rnd(-20, 20)), v = P((db) rnd(-20, 20), (db) rnd(-20, 20));
            p2 w = P((db) rnd(-20, 20), (db) rnd(-20, 20)), z = P((db) rnd(-20, 20), (db) rnd(-20, 20));
            seg s1 = {u, v}, s2 = {w, z};
            if(dd(u, v) < 1 || dd(w, z) < 1 || !chkll(s1, s2)) continue;
            p2 p = isll(s1, s2);
            db sc1 = 1 + dis(v - u), sc2 = 1 + dis(z - w);
            if(!eqd(cross(s1, p), 0, 1e-9L * sc1 * (1 + dis(p - u)))) pb.hit("isll 不在第一条直线上", cross(s1, p), 0);
            if(!eqd(cross(s2, p), 0, 1e-9L * sc2 * (1 + dis(p - w)))) pb.hit("isll 不在第二条直线上", cross(s2, p), 0);
            if(chkss_s(s1, s2) && (!ons(s1, p) || !ons(s2, p))) pb.hit("规范相交时交点不在两段上", ons(s1, p), 1);
        }
        CHECK(pb.bad == 0, "isll 交点同时在两条直线上,规范相交时落在两段内(2 万组)");
    }
    {
        // 退化:退化成点的线段被 chkss 当作「与任何有投影重叠的线段相交」
        printf("  [note] 退化成点的线段上 chkss 只看 crossop == 0(eps = 1e-10 绝对容差)与投影包围盒(不计入失败):"
               "chkss({(2,1e-11),(2,1e-11)}, {(0,0),(4,0)}) 实测 %d —— 该点在 x 轴上方 1e-11,并不是线段上的点,"
               "却因 crossop 归零而被判相交(把 y 抬到 4e-11 以上才变成 %d)。\n",
               (int) chkss(seg{P(2, 1e-11L), P(2, 1e-11L)}, seg{P(0, 0), P(4, 0)}),
               (int) chkss(seg{P(2, 4e-11L), P(2, 4e-11L)}, seg{P(0, 0), P(4, 0)}));
    }

    // ===== 8. disss(返回平方距离) =====
    CHECK(eqd(disss({P(0, 0), P(1, 0)}, {P(0, 1), P(1, 1)}), 1), "disss 平行线段间隔为 1");
    CHECK(eqd(disss({P(0, 0), P(2, 2)}, {P(0, 2), P(2, 0)}), 0), "disss 相交为 0");
    CHECK(eqd(disss({P(0, 0), P(1, 0)}, {P(3, 2), P(4, 2)}), 8), "disss 端点间最短(2^2+2^2 = 8)");
    CHECK(eqd(disss({P(0, 0), P(1, 0)}, {P(3, 4), P(4, 4)}), 20), "disss 端点间最短(2^2+4^2 = 20)");
    CHECK(eqd(disss({P(0, 0), P(4, 0)}, {P(2, -1), P(2, 1)}), 0), "disss 相交为 0(垂直线段)");
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
    {
        Probe pb;  // 距离的性质:对称、非负、相交为 0、被任意端点距离「上界」约束
        For(t, 1, 20000) {
            p2 u = P((db) rnd(-30, 30), (db) rnd(-30, 30)), v = P((db) rnd(-30, 30), (db) rnd(-30, 30));
            p2 w = P((db) rnd(-30, 30), (db) rnd(-30, 30)), z = P((db) rnd(-30, 30), (db) rnd(-30, 30));
            if(dd(u, v) < 1 || dd(w, z) < 1) continue;
            db d = disss(seg{u, v}, seg{w, z});
            if(d < -1e-9L) pb.hit("disss 为负", d, 0);
            if(!eqd(d, disss(seg{w, z}, seg{u, v}), 1e-9L * (1 + d))) pb.hit("disss 不对称", d, disss(seg{w, z}, seg{u, v}));
            db ub = min(min(dd(u, w), dd(u, z)), min(dd(v, w), dd(v, z)));
            if(d > ub + 1e-6L * (1 + ub)) pb.hit("disss 大于端点最小距离", d, ub);
            if(chkss(seg{u, v}, seg{w, z}) && d != 0) pb.hit("相交却距离非 0", d, 0);
            if(refChkss(u, v, w, z) && d > 1e-6L * (1 + ub)) pb.hit("规范相交却距离非 0", d, 0);
        }
        CHECK(pb.bad == 0, "disss 对称/非负/不超过端点距离/相交为 0(2 万组)");
    }

    // ===== 9. area / contain(0 外 / 1 边界 / 2 内)=====
    p2 sq[4] = {P(0, 0), P(2, 0), P(2, 2), P(0, 2)};
    p2 sqr[4] = {P(0, 0), P(0, 2), P(2, 2), P(2, 0)};
    CHECK(eqd(area(4, sq), 4) && eqd(area(4, sqr), -4), "area 带号面积(逆时针为正、顺时针为负)");
    CHECK(eqd(area(2, sq), 0) && eqd(area(1, sq), 0) && eqd(area(0, sq), 0), "area n<=2 为 0");
    {
        int bad = 0;
        For(x, -2, 6) For(y, -2, 6) {
            int got = contain(4, sq, P((db) x / 2, (db) y / 2));
            bool in = x >= 0 && x <= 4 && y >= 0 && y <= 4;
            bool bd = in && (x == 0 || x == 4 || y == 0 || y == 4);
            int want = bd ? 1 : (in ? 2 : 0);
            if(got != want) ++bad;
        }
        CHECK(bad == 0, "contain 正方形网格 0(外)/1(边界)/2(内)全部正确(9×9 网格,含顶点/边上点/角外)");
    }
    CHECK(contain(4, sqr, P(1, 1)) == 2 && contain(4, sqr, P(5, 5)) == 0 && contain(4, sqr, P(0, 1)) == 1, "contain 对顺时针多边形同样正确(只看几何不看旋向)");
    {
        // 凹多边形(梳子形)+ 边界点:与独立射线法对照,网格逐点扫
        vector<p2> comb = {P(0, 0), P(6, 0), P(6, 1), P(5, 1), P(5, 5), P(4, 5), P(4, 1), P(3, 1), P(3, 5), P(2, 5), P(2, 1), P(1, 1), P(1, 5), P(0, 5)};
        Probe pb;
        For(x, -1, 7) For(y, -1, 6) {
            p2 q = P(x, y);
            int got = contain(comb.size(), comb.data(), q), want = refContain(comb, q);
            if(got != want) {
                pb.hitp("凹多边形 contain 与独立射线法不一致", P(got, 0), P(want, 0));
                if(pb.bad) break;
            }
        }
        CHECK(pb.bad == 0, "contain 凹多边形(梳子形)逐点与独立射线法一致(含凹口、边界、角点,8×7 网格)");
        CHECK(eqd(Abs(area(comb.size(), comb.data())), 18), "area 凹多边形手算面积 = 18(梳子形:底 6×1 + 4 个 1×3 的齿)");
    }
    {
        // 随机星形(关于原点可见)多边形:与独立射线法大量对拍
        Probe pb;
        int nt = 0;
        For(t, 1, 300) {
            int n = (int) rnd(3, 12);
            vector<db> ang(n), rr(n);
            ForD(i, 0, n) ang[i] = 2 * pi * i / n + (db) rnd(-100, 100) / 1000.0;
            ForD(i, 0, n) rr[i] = (db) rnd(100, 1000) / 10;
            vector<p2> poly(n);
            ForD(i, 0, n) poly[i] = P(rr[i] * cos(ang[i]), rr[i] * sin(ang[i]));
            ForD(it, 0, 60) {
                p2 q = P((db) rnd(-1200, 1200) / 10, (db) rnd(-1200, 1200) / 10);
                if(contain(n, poly.data(), q) != refContain(poly, q)) { ++pb.bad; ++nt; }
            }
            if(!eqd(Abs(area(n, poly.data())), Abs(refArea(poly)), 1e-6L * (1 + Abs(refArea(poly))))) ++pb.bad;
        }
        CHECK(pb.bad == 0, "contain/area 与独立射线法/shoelace 一致(300 个随机星形多边形 × 60 点 = 1.8 万次)");
    }
    {
        Probe pb;  // 平移/镜像不变性:点集平移不改变包含关系,镜像后 inside 变 outside
        For(t, 1, 200) {
            int n = (int) rnd(3, 10);
            vector<p2> poly(n);
            ForD(i, 0, n) poly[i] = P((db) rnd(-50, 50), (db) rnd(-50, 50));
            if(Abs(area(n, poly.data())) < 1) continue;
            db dx = (db) rnd(-100, 100) / 7, dy = (db) rnd(-100, 100) / 7;
            vector<p2> mv(n);
            ForD(i, 0, n) mv[i] = poly[i] + P(dx, dy);
            ForD(it, 0, 40) {
                p2 q = P((db) rnd(-60, 60), (db) rnd(-60, 60));
                if(contain(n, poly.data(), q) != contain(n, mv.data(), q + P(dx, dy))) ++pb.bad;
            }
        }
        CHECK(pb.bad == 0, "contain 平移不变性(200 个随机多边形 × 40 点)");
    }
    {
        // 退化观察:自交多边形/退化多边形的约定
        vector<p2> bow = {P(0, 0), P(4, 4), P(4, 0), P(0, 4)};  // 蝴蝶结(自交)
        printf("  [note] contain/area 对自交多边形的约定(不计入失败):蝴蝶结 area = %Lg(两个三角形符号相反,相消为 0),"
               "contain(中心 (2,2) 是交叉点) = %d(边界)。模板只保证简单多边形的语义。\n",
               area(4, bow.data()), contain(4, bow.data(), P(2, 2)));
        vector<p2> deg = {P(0, 0), P(1, 0), P(2, 0)};
        printf("  [note] 全共线「多边形」:contain 恒为 %d(0/2 都不是),area = %Lg。\n", contain(3, deg.data(), P(1, 0)), area(3, deg.data()));
    }

    // ===== 10. convex_hull / convex_diameter =====
    {
        // 10.1 小整数坐标随机点集(与独立 Andrew 参考逐位对照 + 全部几何性质)
        Probe pb;
        int dupcnt = 0;
        For(t, 1, 3000) {
            int n = (int) rnd(1, 14), R = (int) rnd(1, 8);
            vector<p2> v;
            if(t <= 1500) {  // 前半允许大量重复点(每个坐标只有 (2R+1)^2 种可能)
                ForD(i, 0, n) v.push_back(P((db) rnd(-R, R), (db) rnd(-R, R)));
            } else {  // 后半保证互不相同,覆盖一般位置
                ForD(i, 0, n) v.push_back(P((db) rnd(-60, 60), (db) rnd(-60, 60)));
                ForD(i, 0, n) ForD(j, 0, (int) v.size()) if(i != j && refCmpP(v[i], v[j]) == 0) { v[i].x += 1000; }
            }
            if(!hullCheckOne(v, true)) {
                if(!pb.bad) {  // 第一条反例:把输入点集打出来
                    printf("      反例输入(%d 个点):", n);
                    ForD(i, 0, n) pr(v[i]);
                    printf("\n");
                }
                ++pb.bad, ++dupcnt;
            }
        }
        CHECK(pb.bad == 0, "convex_hull 与独立 Andrew 参考一致(3000 组:随机网格 + 互异点;含 k<=n、逆时针、"
                           "凸性、首点为最小点、包含全部输入点、顶点序列逐一相同)");
    }
    {
        // 10.2 输入点数的上界与「顶点必须来自输入点」
        Probe pb;
        For(t, 1, 3000) {
            int n = (int) rnd(1, 40);
            vector<p2> v;
            ForD(i, 0, n) v.push_back(P((db) rnd(-1000000, 1000000), (db) rnd(-1000000, 1000000)));
            p2 a2[64], h[64];
            ForD(i, 0, n) a2[i] = v[i];
            int k = convex_hull(n, a2, h);
            if(k > n) { pb.hit("顶点数超过输入点数", k, n); continue; }
            if(n == 1) {  // n==1 时模板不写 b(第 10 段 [note]),只查返回值
                if(k != 1) pb.hit("n=1 返回值应为 1", k, 1);
                continue;
            }
            ForD(i, 0, k) {  // 排序会重排 a2,所以和 a2 比(而非原始 v)
                bool found = false;
                ForD(j, 0, n) if(refCmpP(h[i], a2[j]) == 0) found = true;
                if(!found) pb.hit("凸包顶点不是输入点", 0, 0);
            }
            // 凸包面积 <= 任何内接多边形面积不可能更大:用 shoelace 面积非负(逆时针)
            if(k >= 3 && refArea(vector<p2>(h, h + k)) < -1e-6L * (1 + Abs(refArea(vector<p2>(h, h + k))))) pb.hit("凸包面积(带号)为负", 0, 0);
        }
        CHECK(pb.bad == 0, "convex_hull 性质:k <= n、顶点全都来自输入点、带号面积非负(3000 组,坐标 ±1e6)");
    }
    {
        // 10.3 全部共线 / 全部同点 / n = 0,1,2,3 / 圆上点 / 大量共线点
        vector<p2> line;
        For(i, -50, 50) line.push_back(P(i, 2 * i + 3));
        CHECK(hullCheckOne(line), "凸包:101 个共线点 -> 只有两个端点");
        {
            vector<p2> allsame(37, P(-8, 11));
            p2 a2[64], h[64];
            ForD(i, 0, 37) a2[i] = allsame[i];
            int k = convex_hull(37, a2, h);
            CHECK(k == 2 && eqp(h[0], P(-8, 11)) && eqp(h[1], P(-8, 11)), "凸包:37 个完全相同的点 -> 返回 2 个重复的同一个点(k=2,不是 1)");
        }
        {
            p2 one[1] = {P(7, 8)}, h[8] = {P(-1, -1), P(-1, -1), P(-1, -1), P(-1, -1)};
            int k = convex_hull(1, one, h);
            printf("  [note] convex_hull(n<=1) 返回 %d 但不写 b(n=1 时 b[0] 仍是哨兵 (-1,-1)=%d;n=0 也返回 %d)"
                   " —— 调用方必须在 n<=1 时自己处理(不计入失败)。\n",
                   k, (int) eqp(h[0], P(-1, -1)), convex_hull(0, one, h));
        }
        {
            vector<p2> two = {P(3, 4), P(0, 0)};
            CHECK(hullCheckOne(two), "凸包:n=2 返回两个端点(顺序为字典序)");
            vector<p2> tri = {P(0, 0), P(4, 0), P(2, 3)};
            CHECK(hullCheckOne(tri), "凸包:n=3 非共线 -> 三个顶点,逆时针");
            vector<p2> tri2 = {P(0, 0), P(1, 1), P(2, 2)};
            CHECK(hullCheckOne(tri2), "凸包:n=3 共线 -> 两个端点");
            vector<p2> rep = {P(1, 1), P(1, 1), P(1, 1), P(1, 1), P(1, 1)};
            CHECK(hullCheckOne(rep), "凸包:5 个同点 -> 两个重复点");
            vector<p2> mixed = {P(0, 0), P(0, 0), P(3, 4), P(3, 4), P(1, 2)};
            CHECK(hullCheckOne(mixed), "凸包:带重复点的共线点集 -> 两个端点");
        }
        {
            // 圆上 40 个点(精度受 cos/sin 支配,只做性质级断言)
            vector<p2> cir;
            ForD(i, 0, 40) {
                db ang = 2 * pi * i / 40;
                cir.push_back(P(1000 * cos(ang), 1000 * sin(ang)));
            }
            p2 a2[64], h[64];
            ForD(i, 0, 40) a2[i] = cir[i];
            int k = convex_hull(40, a2, h);
            Probe pb;
            if(k != 40) pb.hit("圆上 40 点凸包顶点数应为 40", k, 40);
            ForD(i, 0, k) {
                p2 u = h[i], w = h[(i + 1) % k];
                if(refDet(w - u, cir[(2 * i) % 40] - u) < -1e-6L) pb.hit("圆上点凸包非凸/漏点", 0, 0);
            }
            db ar = Abs(refArea(vector<p2>(h, h + k)));
            if(!eqd(ar, pi * 1e6L * sin(2 * pi / 40) * 40 / (2 * pi), 1e-6 * (1 + ar))) pb.hit("圆内接正 40 边形面积不对", ar, pi * 1e6L * sin(2 * pi / 40) * 40 / (2 * pi));
            CHECK(pb.bad == 0, "凸包:半径 1000 圆上 40 点 -> 全部 40 个顶点,逆时针凸,面积 = (1/2)·40·R²·sin(2π/40)");
        }
        {
            // 大量共线点 + 少量角点:去共线语义必须只留角点
            Probe pb;
            For(t, 1, 400) {
                int kk = (int) rnd(0, 12);
                vector<p2> v;
                ForD(i, 0, kk) v.push_back(P(0, i));  // 左边上的共线点
                ForD(i, 0, kk) v.push_back(P(100, i));  // 右边上的共线点
                ForD(i, 0, kk) v.push_back(P(i, 0));
                v.push_back(P(0, 0)), v.push_back(P(100, 0)), v.push_back(P(100, 100)), v.push_back(P(0, 100));
                if(!hullCheckOne(v)) ++pb.bad;
                p2 a2[64], h[64];
                ForD(i, 0, (int) v.size()) a2[i] = v[i];
                int got = convex_hull(v.size(), a2, h);
                if(got != 4) pb.hit("共线点应被去掉,只剩 4 个角点", got, 4);
            }
            CHECK(pb.bad == 0, "凸包 nos=0:矩形边上塞满共线点 -> 仍只剩 4 个角点(400 组,每边 0~12 个点)");
        }
        {
            vector<p2> pts = {P(0, 0), P(1, 0), P(2, 0), P(3, 0), P(3, 3), P(0, 3), P(0, 1), P(1, 1)};
            vector<p2> t0 = pts, t1 = pts;
            p2 h0[32], h1[32];
            int k0 = convex_hull(t0.size(), t0.data(), h0, 0), k1 = convex_hull(t1.size(), t1.data(), h1, 1);
            CHECK(k0 == 4 && eqp(h0[0], P(0, 0)) && eqp(h0[1], P(3, 0)) && eqp(h0[2], P(3, 3)) && eqp(h0[3], P(0, 3)), "convex_hull nos=0 去掉共线点(只剩 4 个角点)");
            CHECK(k1 == 7, "convex_hull nos=1 保留共线点(4 角点 + 3 个边上点)");
        }
        {
            // nos=1 的语义:凸包边界上的输入点一个不少
            Probe pb;
            For(t, 1, 300) {
                int n = (int) rnd(4, 12);
                vector<p2> v;
                ForD(i, 0, n) v.push_back(P((db) rnd(-6, 6), (db) rnd(-6, 6)));
                p2 a2[64], h[64];
                ForD(i, 0, n) a2[i] = v[i];
                int k = convex_hull(n, a2, h, 1);
                vector<p2> H(h, h + k);
                ForD(i, 0, n) {  // 每个输入点都必须仍然「在凸包内」(容差级);nos=1 只是不去掉共线点
                    bool in = k < 3;
                    ForD(j, 0, k) {
                        db c = refDet(h[(j + 1) % k] - h[j], v[i] - h[j]);
                        if(c >= -1e-9L * (1 + Abs(v[i].x) + Abs(v[i].y))) in = true;
                    }
                    if(!in) {
                        if(!pb.bad) {
                            printf("      反例:输入点 "), pr(v[i]), printf(" 在 nos=1 凸包外侧;hull(%d 点):", k);
                            ForD(j, 0, k) pr(h[j]);
                            printf("\n");
                        }
                        pb.hit("nos=1 时凸包不再包含输入点", 0, 0);
                    }
                }
                if(k > 2 * n) pb.hit("nos=1 顶点数异常膨胀(k > 2n)", k, n);  // 全共线时模板会来回走一遍,k 可达 2n-2
            }
            CHECK(pb.bad == 0, "convex_hull nos=1:仍然包含所有输入点、顶点数不膨胀(300 组)");
            printf("  [note] nos=1 的「保留共线点」判据是 crossop == 0,门槛是 eps = 1e-10 的**绝对**叉积;"
                   "long double 下整数坐标(<= 1e9)的整数叉积是精确的,所以这一档没有踩到(不计入失败),"
                   "但一旦坐标本身带舍入(如来自旋转/交点的浮点结果),绝对 eps 就会把「离边不到 1e-10 的点」也算作在边上。\n");
        }
    }
    {
        // 10.4 旋转不变性 / 平移不变性 / 缩放不变性(几何性质,不逐位比较)
        Probe pb;
        For(t, 1, 500) {
            int n = (int) rnd(3, 12);
            vector<p2> v;
            for(int i = 0; i < n; i++) v.push_back(P((db) rnd(-100, 100), (db) rnd(-100, 100)));
            vector<p2> R = refHull(v);
            p2 a2[64], h[64];
            ForD(i, 0, n) a2[i] = v[i];
            int k = convex_hull(n, a2, h);
            if(k != (int) R.size()) { pb.hit("原始点数与参考不符", k, R.size()); continue; }
            db s0 = Abs(refArea(R));
            // 旋转
            db ang = (db) rnd(-3000, 3000) / 1000;
            vector<p2> rv(n), rH;
            ForD(i, 0, n) rv[i] = rot(v[i], ang);
            ForD(i, 0, n) a2[i] = rv[i];
            int k2 = convex_hull(n, a2, h);
            if(k2 != (int) R.size()) pb.hit("旋转后顶点数变了", k2, R.size());
            db s1 = Abs(refArea(vector<p2>(h, h + k2)));
            if(!eqd(s1, s0, 1e-6 * (1 + s0))) pb.hit("旋转后面积变了", s1, s0);
            // 平移
            ForD(i, 0, n) a2[i] = v[i] + P(1234.5L, -987.25L);
            int k3 = convex_hull(n, a2, h);
            if(k3 != (int) R.size()) pb.hit("平移后顶点数变了", k3, R.size());
            if(!eqd(Abs(refArea(vector<p2>(h, h + k3))), s0, 1e-6 * (1 + s0))) pb.hit("平移后面积变了", 0, 0);
            // 正数缩放
            ForD(i, 0, n) a2[i] = v[i] * 7.0L;
            int k4 = convex_hull(n, a2, h);
            if(k4 != (int) R.size()) pb.hit("缩放后顶点数变了", k4, R.size());
            if(!eqd(Abs(refArea(vector<p2>(h, h + k4))), s0 * 49, 1e-5 * (1 + s0 * 49))) pb.hit("缩放后面积不是 49 倍", Abs(refArea(vector<p2>(h, h + k4))), s0 * 49);
        }
        CHECK(pb.bad == 0, "convex_hull 不变量:旋转/平移保面积与顶点数,7 倍缩放面积 ×49(500 组)");
    }
    {
        // 10.5 convex_diameter:暴搜最远点对(返回平方距离)
        Probe pb;
        For(t, 1, 3000) {
            int n = (int) rnd(2, 12);
            vector<p2> v;
            ForD(i, 0, n) v.push_back(P((db) rnd(-70, 70), (db) rnd(-70, 70)));
            p2 a2[64], h[64];
            ForD(i, 0, n) a2[i] = v[i];
            int k = convex_hull(n, a2, h);
            db bf = 0;
            ForD(i, 0, k) ForD(j, 0, k) bf = max(bf, dd(h[i], h[j]));
            db got = convex_diameter(k, h);
            if(!eqd(got, bf, 1e-9 * (1 + bf))) pb.hit("convex_diameter != 暴搜最远点对", got, bf);
            if(!eqd(dis2(a2[0] - a2[1]), 0, 1e300L)) {}  // no-op 保序
        }
        p2 two[2] = {P(0, 0), P(3, 4)};
        if(!eqd(convex_diameter(2, two), 25)) pb.hit("n=2 直径", convex_diameter(2, two), 25);
        p2 col[3] = {P(0, 0), P(1, 1), P(3, 3)};
        if(!eqd(convex_diameter(3, col), 18)) pb.hit("全共线直径", convex_diameter(3, col), 18);
        p2 dup[2] = {P(1, 1), P(1, 1)};
        if(!eqd(convex_diameter(2, dup), 0)) pb.hit("两点重合直径", convex_diameter(2, dup), 0);
        p2 one[1] = {P(5, 5)};
        if(!eqd(convex_diameter(1, one), 0)) pb.hit("n=1 直径", convex_diameter(1, one), 0);
        vector<p2> cird;
        ForD(i, 0, 24) cird.push_back(P(1000 * cos(2 * pi * i / 24), 1000 * sin(2 * pi * i / 24)));
        db bf = 0;
        ForD(i, 0, 24) ForD(j, 0, 24) bf = max(bf, dd(cird[i], cird[j]));
        if(!eqd(convex_diameter(24, cird.data()), bf, 1e-9 * (1 + bf))) pb.hit("圆上点直径", convex_diameter(24, cird.data()), bf);
        CHECK(pb.bad == 0, "convex_diameter = 暴搜最远点对(3000 组随机凸包 + n=1/2/全共线/重合点/圆上点的边界用例)");
    }

    // ===== 11. convex_cut =====
    {
        // 11.1 与独立参考半平面裁剪对照(整数/简单浮点)
        Probe pb;
        For(t, 1, 3000) {
            int n = (int) rnd(3, 10);
            vector<p2> v;
            ForD(i, 0, n) v.push_back(P((db) rnd(-30, 30), (db) rnd(-30, 30)));
            p2 a2[64], h[64];
            ForD(i, 0, n) a2[i] = v[i];
            int k = convex_hull(n, a2, h);
            if(k < 3) continue;
            vector<p2> H(h, h + k);
            seg q = {P((db) rnd(-60, 60), (db) rnd(-60, 60)), P((db) rnd(-60, 60), (db) rnd(-60, 60))};
            if(dd(q.x, q.y) < 1) continue;
            vector<p2> got = convex_cut(k, H.data(), q), want = refCut(H, q);
            if(got.size() != want.size()) { pb.hit("裁剪后顶点数与参考不同", got.size(), want.size()); continue; }
            ForD(i, 0, want.size()) if(!eqp(got[i], want[i], 1e-6)) pb.hitp("裁剪顶点与参考不同", got[i], want[i]);
            if(!eqd(cutArea(got), cutArea(want), 1e-6 * (1 + cutArea(want)))) pb.hit("裁剪后面积与参考不同", cutArea(got), cutArea(want));
        }
        CHECK(pb.bad == 0, "convex_cut 与独立参考半平面裁剪逐顶点一致(3000 组随机凸包)");
    }
    {
        // 11.2 性质级:顶点数 <= n、顶点来自输入或落在切割线上、旋向不变、凸性保持、面积不增
        Probe pb;
        For(t, 1, 3000) {
            int n = (int) rnd(3, 12);
            vector<p2> v;
            ForD(i, 0, n) v.push_back(P((db) rnd(-40, 40), (db) rnd(-40, 40)));
            p2 a2[64], h[64];
            ForD(i, 0, n) a2[i] = v[i];
            int k = convex_hull(n, a2, h);
            if(k < 3) continue;
            vector<p2> H(h, h + k);
            seg q = {P((db) rnd(-70, 70), (db) rnd(-70, 70)), P((db) rnd(-70, 70), (db) rnd(-70, 70))};
            if(dd(q.x, q.y) < 1) continue;
            vector<p2> got = convex_cut(k, H.data(), q);
            // 凸 k 边形被一条直线切开,inside 顶点数 i 与跨越边数 c 满足 i + c <= k + 1,故输出点数 <= k + 1;
            // 这里再放宽 1 以吸收交点浮点舍入带来的等价重复点(实测常见退化是 k+1,见本段末尾 [note])
            if((int) got.size() > k + 2) {
                if(!pb.bad) {
                    printf("      反例:输入 %d 边形:", k);
                    ForD(j, 0, k) pr(H[j]);
                    printf("\n            切割线 "), pr(q.x), printf("-"), pr(q.y), printf(" -> 输出 %zu 个点:", got.size());
                    ForD(j, 0, got.size()) pr(got[j]);
                    printf("\n            crossop(q,H[j]) =");
                    ForD(j, 0, k) printf(" %d", crossop(q, H[j]));
                    printf("\n");
                }
                pb.hit("裁剪后顶点数异常", got.size(), k + 1);
            }
            db a0 = Abs(refArea(H)), a1 = cutArea(got);
            if(a1 > a0 + 1e-6 * (1 + a0)) pb.hit("裁剪后面积变大", a1, a0);
            ForD(i, 0, got.size()) {
                if(cross(q, got[i]) < -1e-9L * (1 + dis(got[i] - q.x)) * dis(q.dir())) pb.hit("裁剪后顶点跑到半平面外侧", cross(q, got[i]), 0);
                bool fromIn = false;
                ForD(j, 0, k) if(eqp(got[i], H[j], 1e-7)) fromIn = true;
                if(!fromIn && Abs(cross(q, got[i])) > 1e-6L * (1 + dis(got[i] - q.x)) * dis(q.dir())) pb.hit("裁剪顶点既不是输入顶点也不在切割线上", 0, 0);
            }
            ForD(i, 0, got.size()) {
                p2 u = got[i], w = got[(i + 1) % got.size()], z = got[(i + 2) % got.size()];
                if(got.size() >= 3 && refDet(w - u, z - u) < -1e-6L) pb.hit("裁剪后非凸", 0, 0);
            }
            if(got.size() >= 3 && refArea(got) < -1e-7L * (1 + a0)) pb.hit("裁剪后旋向反转", refArea(got), 0);
        }
        CHECK(pb.bad == 0, "convex_cut 性质:顶点数 <= n+2、面积不增、顶点只能来自输入顶点或切割线、保持凸性与逆时针(3000 组)");
        printf("  [note] convex_cut 的输出点数可能达到 n+1(不计入失败):k 边形与直线相交时 inside 顶点数 i 与跨越边数 c 满足\n"
               "         i + c = k + 1 或 k,所以 k+1 是正常上界;并且当某条边的两个端点 crossop 都为 1 却仍被 isll 判出交点时,\n"
               "         a1+a2 相消会留下 1e-16 量级的偏差,使交点落在半平面外一丝丝(实测 cross = -1.1e-16,crossop 仍为 0)。\n");
        printf("  [note] 同一个半平面把切割线端点顺序反过来,输出点数会变(k=4 的实测:5 个点 vs 4 个点),面积差 < 1e-14 —— "
               "这是等价点未合并造成的,不是逻辑错误,但按「点数」比较时要当心。\n");
    }
    {
        // 11.3 解析用例
        p2 sq4[4] = {P(0, 0), P(4, 0), P(4, 4), P(0, 4)};
        Probe pb;
        auto S = [&](vector<p2> r) { return cutArea(r); };
        vector<p2> half = convex_cut(4, sq4, seg{P(2, -1), P(2, 5)});
        if(!eqd(S(half), 8)) pb.hit("正方形切一半面积应为 8", S(half), 8);
        vector<p2> tri = convex_cut(4, sq4, seg{P(0, 0), P(4, 4)});
        if(!eqd(S(tri), 8)) pb.hit("沿对角线切面积应为 8", S(tri), 8);
        if(!eqd(S(convex_cut(4, sq4, seg{P(0, 0), P(0, 4)})), 0)) pb.hit("沿边切面积应为 0", 0, 0);
        vector<p2> out = convex_cut(4, sq4, seg{P(5, -1), P(-1, -1)});  // 保留侧在多边形之外 -> 应返回 0 个点
        if((int) out.size() != 0) pb.hit("切掉全部应返回 0 个点", out.size(), 0);
        vector<p2> all = convex_cut(4, sq4, seg{P(-1, 5), P(-1, -1)});  // 半平面完全包含多边形 -> 面积不变
        if(!eqd(S(all), 16)) pb.hit("半平面完全包含多边形时面积不变", S(all), 16);
        vector<p2> all2 = convex_cut(4, sq4, seg{P(-1, -1), P(5, -1)});  // 另一条支持线,同样完全包含
        if(!eqd(S(all2), 16)) pb.hit("另一方向的支持线面积也不变", S(all2), 16);
        vector<p2> touch = convex_cut(4, sq4, seg{P(-1, 0), P(-1, 0)});  // 退化切割线(零长)
        printf("  [note] convex_cut 的切割线若退化成点(零长),返回值无意义:实测 %zu 个点(不计入失败)。\n", touch.size());
        CHECK(pb.bad == 0, "convex_cut 解析用例:切半面积 8、对角线切面积 8、沿边切面积 0、切空返回 0、完全包含面积 16");
    }
    {
        // 11.4 连续两次切割 = 与先合并半平面一致(面积单调不增、可交换)
        Probe pb;
        For(t, 1, 800) {
            int n = (int) rnd(3, 10);
            vector<p2> v;
            ForD(i, 0, n) v.push_back(P((db) rnd(-50, 50), (db) rnd(-50, 50)));
            p2 a2[64], h[64];
            ForD(i, 0, n) a2[i] = v[i];
            int k = convex_hull(n, a2, h);
            if(k < 3) continue;
            vector<p2> H(h, h + k);
            seg q1 = {P((db) rnd(-80, 80), (db) rnd(-80, 80)), P((db) rnd(-80, 80), (db) rnd(-80, 80))};
            seg q2 = {P((db) rnd(-80, 80), (db) rnd(-80, 80)), P((db) rnd(-80, 80), (db) rnd(-80, 80))};
            if(dd(q1.x, q1.y) < 1 || dd(q2.x, q2.y) < 1) continue;
            vector<p2> r1 = refCut(H, q1);           // 参考裁剪保证凸性,才能再喂给 convex_cut
            if(r1.size() < 3) continue;
            vector<p2> a = refCut(r1, q2);
            vector<p2> b2 = convex_cut(r1.size(), r1.data(), q2);
            if(!eqd(cutArea(a), cutArea(b2), 1e-5 * (1 + cutArea(a)))) pb.hit("连续裁剪面积不一致", cutArea(b2), cutArea(a));
        }
        CHECK(pb.bad == 0, "convex_cut 连续两次切割的面积与参考一致(800 组,切割线随机)");
    }

    // ===== 12. circumcircle_diameter / circumcenter / 最小圆覆盖 =====
    {
        auto c = circumcenter(P(0, 0), P(2, 0), P(1, 1));
        CHECK(eqd(circumcircle_diameter(P(0, 0), P(2, 0), P(1, 1)), 2) && eqp(P(c[0], c[1]), P(1, 0)), "外接圆:直角三角形斜边为直径、圆心在斜边中点");
    }
    CHECK(circumcircle_diameter(P(0, 0), P(1, 0), P(2, 0)) == -1, "共线三点返回 -1");
    CHECK(circumcircle_diameter(P(0, 0), P(0, 0), P(1, 1)) == -1 && circumcircle_diameter(P(0, 0), P(1, 1), P(1, 1)) == -1, "重合两点(退化)也返回 -1");
    {
        auto c = circumcenter(P(0, 0), P(1, 0), P(0, 1));
        CHECK(eqp(P(c[0], c[1]), P(0.5L, 0.5L)) && eqd(circumcircle_diameter(P(0, 0), P(1, 0), P(0, 1)), sqrtl(2.0L)), "等腰直角三角形:外心 (1/2,1/2),直径 √2");
        db d3 = circumcircle_diameter(P(0, 0), P(1, 0), P(0.5L, sqrtl(3.0L) / 2));  // 单位等边三角形
        CHECK(eqd(d3, 2 / sqrtl(3.0L), 1e-12), "单位等边三角形外接圆直径 = 2/√3");
    }
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
    {
        Probe pb;  // 外心与自写解析式一致;圆心必须在三角形所在平面内(到三边距离 = 半弦长约束)
        For(t, 1, 20000) {
            p2 A = P((db) rnd(-80, 80), (db) rnd(-80, 80)), B = P((db) rnd(-80, 80), (db) rnd(-80, 80)), C = P((db) rnd(-80, 80), (db) rnd(-80, 80));
            db ar2 = Abs((B - A).det(C - A));
            if(ar2 < 1) continue;
            auto cc = circumcenter(A, B, C);
            p2 O = P(cc[0], cc[1]), O2 = refCircumcenter(A, B, C);
            db R = dis(A - B) * dis(B - C) * dis(C - A) / (2 * ar2);
            if(!eqp(O, O2, 1e-7L * (1 + R))) pb.hitp("circumcenter 与独立解析式不一致", O, O2);
            // 三角形顶点都在以 O 为心 R 为半径的圆上 -> 原点到圆心的距离满足幂等式
            if(!eqd(dis(O - A), R, 1e-8L * (1 + R))) pb.hit("O 到 A 的距离不是 R", dis(O - A), R);
            if(!eqd(dis(O - B), R, 1e-8L * (1 + R))) pb.hit("O 到 B 的距离不是 R", dis(O - B), R);
            if(!eqd(dis(O - C), R, 1e-8L * (1 + R))) pb.hit("O 到 C 的距离不是 R", dis(O - C), R);
            if(!eqd(circumcircle_diameter(A, B, C), 2 * R, 1e-8L * (1 + R))) pb.hit("直径 != 2R", circumcircle_diameter(A, B, C), 2 * R);
            // 直角/钝角/锐角:2R 必须不小于最长边(单调性)
            db mx = max(dis(A - B), max(dis(B - C), dis(C - A)));
            if(2 * R < mx - 1e-7L * (1 + mx)) pb.hit("直径小于最长边", 2 * R, mx);
        }
        CHECK(pb.bad == 0, "circumcenter 与独立解析式一致、三点共圆且直径 = 2R ≥ 最长边(2 万组)");
    }
    {
        Probe pb;  // 最小圆覆盖暴力:三角形外接圆必须覆盖所有点(性质级)
        For(t, 1, 400) {
            int n = (int) rnd(3, 7);
            vector<p2> v;
            ForD(i, 0, n) v.push_back(P((db) rnd(-30, 30), (db) rnd(-30, 30)));
            p2 O = bruteMEC(v);
            db r2 = 0;
            ForD(i, 0, n) r2 = max(r2, dd(O, v[i]));
            ForD(i, 0, n) if(dd(O, v[i]) > r2 + 1e-9L) pb.hit("bruteMEC 未覆盖点", 0, 0);
            // 对三角形:O 必须是某两点中点或三点外心;验证模板外心确实给出这个圆的圆心
            ForD(i, 0, n) ForD(j, i + 1, n) ForD(k, j + 1, n) {
                if(sign(refDet(v[j] - v[i], v[k] - v[i])) == 0) continue;
                auto cc = circumcenter(v[i], v[j], v[k]);
                p2 O2 = P(cc[0], cc[1]);
                db R2 = dis(v[i] - v[j]) * dis(v[j] - v[k]) * dis(v[k] - v[i]) / (2 * Abs(refDet(v[j] - v[i], v[k] - v[i])));
                if(!eqd(dis(O2 - v[i]), R2, 1e-7L * (1 + R2))) pb.hit("三点外接圆半径不对", dis(O2 - v[i]), R2);
            }
        }
        CHECK(pb.bad == 0, "最小圆覆盖暴力 vs circumcenter:三点外接圆半径与到三点距离一致(400 组点集,全部三角形枚举)");
    }

    // ===== 13. 数值抗性:±1e9 量级只做性质级断言 =====
    {
        Probe pb;
        For(t, 1, 4000) {  // 1e9 量级随机点集的凸包
            int n = (int) rnd(3, 12);
            vector<p2> v;
            ForD(i, 0, n) v.push_back(P((db) rnd(-1000000000, 1000000000), (db) rnd(-1000000000, 1000000000)));
            p2 a2[64], h[64];
            ForD(i, 0, n) a2[i] = v[i];
            int k = convex_hull(n, a2, h);
            if(k > n || k < 1) { pb.hit("hull 顶点数越界", k, n); continue; }
            ForD(i, 0, n) ForD(j, 0, k) {  // 凸包必须包含所有输入点(容差 1e-9 相对量级)
                p2 u = h[j], w = h[(j + 1) % k];
                db c = refDet(w - u, v[i] - u);
                if(c < -1e-9L * (1 + Abs(u.x) + Abs(u.y) + Abs(w.x) + Abs(w.y))) pb.hit("1e9 量级凸包漏掉输入点", c, 0);
            }
            db ar = refArea(vector<p2>(h, h + k));
            if(ar < -1e-6L * (1 + Abs(ar))) pb.hit("1e9 量级凸包带号面积为负(旋向错)", ar, 0);
            db bf = 0;
            ForD(i, 0, k) ForD(j, 0, k) bf = max(bf, dd(h[i], h[j]));
            db gd = convex_diameter(k, h);
            if(!eqd(gd, bf, 1e-8L * (1 + bf))) pb.hit("1e9 量级 convex_diameter 与暴搜不符", gd, bf);
            if(gd < -1e-9L) pb.hit("1e9 量级直径为负", gd, 0);
        }
        CHECK(pb.bad == 0, "1e9 量级凸包:包含全部输入点(容差 1e-9 相对)、逆时针面积非负、直径 = 暴搜最远点对(4000 组)");
    }
    {
        Probe pb;
        For(t, 1, 6000) {  // 1e9 量级的线段距离/最近点
            p2 u = P((db) rnd(-1000000000, 1000000000), (db) rnd(-1000000000, 1000000000)), v = u;
            do {
                v = P((db) rnd(-1000000000, 1000000000), (db) rnd(-1000000000, 1000000000));
            } while(dd(u, v) < 1);
            p2 q = P((db) rnd(-1000000000, 1000000000), (db) rnd(-1000000000, 1000000000));
            seg s = {u, v};
            db ref = refNearest(u, v, q), got = nearest(s, q);
            if(!eqd(got, ref, 1e-10L * (1 + ref))) pb.hit("1e9 量级 nearest 与参考不符", got, ref);
            if(got < -1e-8L * (1 + ref)) pb.hit("nearest 为负", got, 0);
            p2 h = proj(s, q);
            db resid = Abs((q - h) * (v - u)) / (dis(v - u) * (1 + dis(q - u)));
            if(resid > 1e-9L) pb.hit("1e9 量级 proj 垂线不垂直", resid, 0);
            if(!eqd(dd(h, q) + dd(h, u), dd(q, u), 1e-8L * (1 + dd(q, u)))) pb.hit("1e9 量级投影破坏勾股关系", 0, 0);
            if(isMid(u, h, v)) {  // 垂足落在线段内:它必须是最近点(ons 在这里会受 eps 影响,故改查「距离为 0」这一等价的容差式)
                if(!eqd(refNearest(u, v, h), 0, 1e-6L * (1 + dd(u, v)))) pb.hit("1e9 量级垂足不在线段上", refNearest(u, v, h), 0);
                db d1 = dd(h, q), d2 = min(dd(u, q), dd(v, q));
                if(d1 > d2 * (1 + 1e-9L) + 1e-6L) pb.hit("1e9 量级垂足不是最近点", d1, d2);
            }
            if(!eqp(reflect(s, reflect(s, q)), q, 1e-7L * (1 + dis(q - u)))) pb.hitp("1e9 量级 reflect 不对合", reflect(s, reflect(s, q)), q);
        }
        CHECK(pb.bad == 0, "1e9 量级 nearest/proj/reflect/ons:参考一致、勾股关系、垂足最短、reflect 对合(6000 组)");
    }
    {
        Probe pb;
        For(t, 1, 8000) {  // 1e9 量级的相交判定(带精确参考,量级内交叉积仍是整数、long double 精确)
            p2 u = P((db) rnd(-1000000, 1000000), (db) rnd(-1000000, 1000000)), v = P((db) rnd(-1000000, 1000000), (db) rnd(-1000000, 1000000));
            p2 w = P((db) rnd(-1000000, 1000000), (db) rnd(-1000000, 1000000)), z = P((db) rnd(-1000000, 1000000), (db) rnd(-1000000, 1000000));
            if(dd(u, v) < 1 || dd(w, z) < 1) continue;
            seg s1 = {u, v}, s2 = {w, z};
            if(chkss(s1, s2) != refChkss(u, v, w, z)) pb.hit("1e6 量级 chkss 与参考不符", chkss(s1, s2), refChkss(u, v, w, z));
            db d = disss(s1, s2), rd = refDisss(u, v, w, z);
            if(!eqd(d, rd, 1e-9L * (1 + rd))) pb.hit("1e6 量级 disss 与参考不符", d, rd);
            if(!chkll(s1, s2)) continue;
            p2 p = isll(s1, s2);
            db sc1 = dis(s1.dir()) * (1 + dis(p - u)), sc2 = dis(s2.dir()) * (1 + dis(p - w));
            if(Abs(cross(s1, p)) > 1e-6L * sc1 || Abs(cross(s2, p)) > 1e-6L * sc2) pb.hit("1e6 量级 isll 交点偏离直线", cross(s1, p), 0);
            if(chkss_s(s1, s2)) {  // 规范相交 -> 交点必须落在两段上(用距离判,避开 ons 的绝对 eps)
                db sc = 1e-6L * (1 + dis(v - u) + dis(z - w));
                if(refNearest(u, v, p) > sc || refNearest(w, z, p) > sc) pb.hit("1e6 量级规范相交时交点不在段上", refNearest(u, v, p), 0);
            }
        }
        CHECK(pb.bad == 0, "1e6 量级 chkss/disss/isll/chkll 与精确参考一致(8000 组)");
    }
    {
        // 1e9 量级 contain / area(整数坐标,det 精确)+ 精确的 gcd 面积公式对照
        Probe pb;
        For(t, 1, 2000) {
            int n = (int) rnd(3, 8);
            vector<p2> v;
            ForD(i, 0, n) v.push_back(P((db) rnd(-1000000000, 1000000000), (db) rnd(-1000000000, 1000000000)));
            p2 a2[64], h[64];
            ForD(i, 0, n) a2[i] = v[i];
            int k = convex_hull(n, a2, h);
            if(k < 3) continue;
            vector<p2> H(h, h + k);
            db A = Abs(area(k, h));
            if(A < -1e-6) pb.hit("1e9 量级 area 为负", A, 0);
            if(!eqd(A, Abs(refArea(H)), 1e-9L * (1 + A))) pb.hit("1e9 量级 area 与 shoelace 不符", A, Abs(refArea(H)));
            ForD(it, 0, 25) {
                p2 q = P((db) rnd(-1000000000, 1000000000), (db) rnd(-1000000000, 1000000000));
                int got = contain(k, H.data(), q), want = refContain(H, q);
                if(got != want) pb.hit("1e9 量级 contain 与独立射线法不符", got, want);
            }
            // 凸包内部点(顶点加权平均)必判为「内」
            p2 c = P(0, 0);
            ForD(i, 0, k) c = c + H[i] * (1.0L / k);
            if(contain(k, H.data(), c) == 0 && k >= 3 && A > 1e9) pb.hit("凸包重心被判为外", 0, 2);
        }
        CHECK(pb.bad == 0, "1e9 量级 area/contain 与独立参考一致(2000 组凸包 × 25 点;含重心必为内部的断言)");
    }
    {
        // 1e9 量级 rot:保长/保角(只做性质断言)
        Probe pb;
        For(t, 1, 4000) {
            p2 x = P((db) rnd(-1000000000, 1000000000), (db) rnd(-1000000000, 1000000000));
            db ang = (db) rnd(-1000000, 1000000) / 100000;
            p2 r = rot(x, ang);
            if(!eqd(dis(r), dis(x), 1e-12L * (1 + dis(x)))) pb.hit("1e9 量级 rot 不保长", dis(r), dis(x));
            if(!eqd(remainderl(r.alpha() - x.alpha() - ang, 2 * pi), 0, 1e-9)) pb.hit("1e9 量级 rot 角度不对", r.alpha() - x.alpha(), ang);
            if(!eqd(x * r, dis2(x) * cos(ang), 1e-9L * dis2(x))) pb.hit("1e9 量级 rot 点积公式不对", x * r, dis2(x) * cos(ang));
            if(!eqp(rot(rot(x, ang), -ang), x, 1e-6L * (1 + dis(x)))) pb.hitp("1e9 量级 rot 逆旋转不还原", rot(rot(x, ang), -ang), x);
        }
        CHECK(pb.bad == 0, "1e9 量级 rot:保长、角度加成、与 cos 的点积公式、逆旋转还原(4000 组)");
    }

    // ===== 14. 大坐标下的圆相关 =====
    {
        Probe pb;
        For(t, 1, 6000) {
            p2 A = P((db) rnd(-1000000000, 1000000000), (db) rnd(-1000000000, 1000000000));
            p2 B = P((db) rnd(-1000000000, 1000000000), (db) rnd(-1000000000, 1000000000));
            p2 C = P((db) rnd(-1000000000, 1000000000), (db) rnd(-1000000000, 1000000000));
            db d = Abs((B - A).det(C - A));
            if(d < 1e6L) continue;  // 过滤近共线的病态三角形(模板 / 参考都无法保证)
            db D = circumcircle_diameter(A, B, C);
            if(D < 0) { pb.hit("非共线却不返回正直径", D, 1); continue; }
            auto cc = circumcenter(A, B, C);
            p2 O = P(cc[0], cc[1]);
            db R = D / 2;
            if(!eqd(dis(O - A), R, 1e-7L * (1 + R)) || !eqd(dis(O - B), R, 1e-7L * (1 + R)) || !eqd(dis(O - C), R, 1e-7L * (1 + R)))
                pb.hit("大坐标外心到三点距离不是 R", dis(O - A), R);
            db mx = max(dis(A - B), max(dis(B - C), dis(C - A)));
            if(D < mx * (1 - 1e-9L)) pb.hit("大坐标直径小于最长边", D, mx);
            // D = abc/(2|det|):与用三角形面积的解析式对照(相对容差)
            db D2 = 2 * dis(A - B) * dis(B - C) * dis(C - A) / (2 * d);
            if(!eqd(D, D2, 1e-9L * (1 + D2))) pb.hit("大坐标直径与解析式不符", D, D2);
            // 外接圆必须覆盖三角形三点(即半径 R 全等)
            if(dis(O - A) > R * (1 + 1e-9L) || dis(O - B) > R * (1 + 1e-9L) || dis(O - C) > R * (1 + 1e-9L))
                pb.hit("大坐标外接圆未覆盖三点", 0, 0);
        }
        CHECK(pb.bad == 0, "1e9 量级外接圆:直径 = abc/(2|det|)、圆心到三点等距为 R、直径 ≥ 最长边(6000 组)");
    }
    {
        Probe pb;
        For(t, 1, 500) {  // 大坐标 + 圆上点凸包
            int n = (int) rnd(8, 16);
            p2 O = P((db) rnd(-1000000000, 1000000000), (db) rnd(-1000000000, 1000000000));
            db R = (db) rnd(1000000, 1000000000);
            vector<p2> v;
            ForD(i, 0, n) {
                db ang = 2 * pi * i / n + (db) rnd(-1000, 1000) / 1000000.0;
                v.push_back(O + P(R * cos(ang), R * sin(ang)));
            }
            p2 a2[64], h[64];
            ForD(i, 0, n) a2[i] = v[i];
            int k = convex_hull(n, a2, h);
            if(k != n) pb.hit("圆上点凸包漏顶点", k, n);
            ForD(i, 0, n) ForD(j, 0, k) {
                p2 u = h[j], w = h[(j + 1) % k];
                if(refDet(w - u, v[i] - u) < -1e-9L * (R * (1 + Abs(R)))) pb.hit("圆上点被判在凸包外", 0, 0);
            }
            db bf = 0;
            ForD(i, 0, k) ForD(j, 0, k) bf = max(bf, dd(h[i], h[j]));
            if(!eqd(convex_diameter(k, h), bf, 1e-8L * (1 + bf))) pb.hit("大坐标圆上点直径不对", convex_diameter(k, h), bf);
        }
        CHECK(pb.bad == 0, "1e9 量级 + 半径 1e6~1e9 圆上点:凸包不漏点、直径 = 暴搜(500 组)");
    }

    // ===== 15. 边界观察:convex_hull 对「全部点相同」的退化输入(n<=1 不写 b 见上)=
    {
        printf("  [note] convex_hull 对完全重合的点集返回 k=2 且两个顶点相同(k=2 而非 1):37 个相同点 -> %d;"
               "但模板对「含重复点的其它点集」是正确的(如 (0,0)×2+(2,0)+(2,2) -> 3 个顶点)。调用方需自行去重。\n",
               [] {
                   p2 s[3] = {P(5, 5), P(5, 5), P(5, 5)}, o[8];
                   return convex_hull(3, s, o);
               }());
    }

    // ===== 16. operator+= (geo.cpp:32,正确) =====
    {
        int bad = 0;
        For(t, 1, 20000) {
            p2 x = P((db) rnd(-1000, 1000), (db) rnd(-1000, 1000)), y = P((db) rnd(-1000, 1000), (db) rnd(-1000, 1000));
            p2 x0 = x, r = (x += y);
            if(!eqp(x, x0 + y) || !eqp(r, x0 + y)) ++bad;
        }
        CHECK(bad == 0, "operator+=(p2&,p2) 语义正确:a += b 后 a == a + b(geo.cpp:32,2 万组随机)");
    }

    // ===== 17. 复合赋值 -= 与 /=:当前 geo.cpp:33-34 把这两个都写成了加法 =====
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
        // 实际后果:任何用 `rel -= ...` 取差向量再算叉积的几何代码,算出来的叉积符号是错的
        vector<p2> gset = {P(0, 0), P(3, 0), P(4, 2), P(3, 4), P(0, 4), P(-1, 2), P(1, 1), P(2, 3)};
        p2 gp[16], gb[32];
        ForD(i, 0, (int) gset.size()) gp[i] = gset[i];
        int gk = lowerChainWithTemplateMinusEq(gset.size(), gp, gb);
        printf("        实际后果(用模板的 -= 求相邻点差向量再 det):%s(正确凸包有 %zu 个顶点、面积 %Lg)\n",
               gk < 0 ? "diff = b[k-1] -= b[k-2] 得到的叉积与真实值不同 -> 整个扫描判据失效" : "叉积与真实值一致",
               refHull(gset).size(), Abs(refArea(refHull(gset))));
        if(!minus_ok) {
            printf("        ==> a -= b 得到的是 a + b,任何用 -= 做差向量的几何代码(凸包/扫描线/极角排序)都会算错\n");
        }
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
        if(!diveq_ok) {
            printf("        ==> 模板里 unit() 用的是 x / dis(x) 没问题,但任何写 a /= d 的代码要么编译不过、要么结果变成加法\n");
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
