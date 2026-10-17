// 最近点对 自测:独立hypot全点对参考、全指数/平方下溢/大平移3000组、随机/退化/20万点
//
// 对拍方式(文件内自写,不用模板):
//   bruteClosest(v):O(n²) 双重循环，用标准 hypot 独立计算距离，不调用模板 dis。
//   另外两条性质级断言:① 答案必须 <= 任意一对具体点的距离;② 答案 == 0 当且仅当存在重复点。
//
// 用例规模:
//   1 解析用例:n=1(无解)、n=2、n=3、重复点、全共线、正方形网格、±1e9 端点
//   2 4000 组 n=2..40 小整数坐标(±8,大量重复点,含 n 很小与大量重合)-> 与暴力逐位对拍
//   3 3000 组 n=2..150 随机坐标(±1e6,含圆上点/细长分布)-> 与暴力对拍(相对 1e-9)
//   4 性质级 3000 组:答案 <= 任意点对距离、答案 == 0 <=> 有重复点、答案非负
//   5 大规模:
//        n = 200000 的 1e4 间隔晶格 + 一对距离 3.5 的点(真答案已知,且这一对横跨中位线,
//        专门测「跨分界线的窄条」逻辑) -> 必须等于 3.5
//        n = 200000 全同点 -> 必须等于 0(并且要跑得快,证明不是 O(n²))
//   6 ±1e9 大坐标 2000 组:与暴力对拍(相对容差)+ 性质断言
#include "../../_check_base.hpp"
#include "二维向量.cpp"
#include "最近点对.cpp"

// 构造点。
static p2 P(db x, db y) { return {x, y}; }
// 绝对值。
static db Abs(db x) { return x < 0 ? -x : x; }
// 带容差的标量判等。
static bool eqd(db a, db b, db tol) { return Abs(a - b) < tol; }
// 打印点(诊断用)。
static void pr(p2 a) { printf("(%g,%g)", (double) a.x, (double) a.y); }
// 打印点集(诊断用,最多前 20 个)。
static void prv(const char *s, const vector<p2> &v) {
    printf("      %s(%zu):", s, v.size());
    for(int i = 0; i < (int) min<size_t>(v.size(), 20); ++i) printf(" "), pr(v[i]);
    printf("%s\n", v.size() > 20 ? " ..." : "");
}

// ---------- 独立参考实现:O(n²) 暴力 ----------
// O(n^2) 暴力最近距离;n < 2 返回 -1。
static db bruteClosest(const vector<p2> &v) {
    int n = v.size();
    if(n < 2) return -1;
    db ans = 1e100L;
    for(int i = 0; i < n; ++i) for(int j = i + 1; j < n; ++j)
        ans = min(ans, hypot(v[i].x - v[j].x, v[i].y - v[j].y));
    return ans;
}

static bool exponent_ranges() {
    mt19937_64 gen(2026101703);
    vector<int> exponents{numeric_limits<db>::min_exponent - numeric_limits<db>::digits + 3,
                          numeric_limits<db>::min_exponent * 3 / 4, -100, 0, 100, 300};
    for(int d = -4; d <= 4; ++d) {
        exponents.push_back(numeric_limits<db>::min_exponent / 2 + d);
        exponents.push_back((numeric_limits<db>::min_exponent - numeric_limits<db>::digits) / 2 + d);
    }
    for(int t = 0; t < 3000; ++t) {
        int e = exponents[t % exponents.size()], n = 2 + gen() % 30;
        db scale = ldexp(db(1), e);
        db shift = t % 3 == 0 ? ldexp(scale, numeric_limits<db>::digits - 8) : 0;
        vector<p2> v;
        for(int i = 0; i < n; ++i)
            v.push_back({shift + db(int(gen() % 201) - 100) * scale,
                         -shift + db(int(gen() % 201) - 100) * scale});
        if(t % 7 == 0) v.push_back(v[0]);
        // 极小点对与相距极远的点混用，不能因远处点的平方溢出而丢掉局部解。
        if(t % 5 == 0) v.push_back({1e90L, -1e90L});
        db got = dis_closest(v), want = bruteClosest(v);
        db tolerance = 16 * numeric_limits<db>::epsilon() * want + 4 * numeric_limits<db>::denorm_min();
        if(abs(got - want) > tolerance || (got == 0) != (want == 0)) {
            printf("closest exponent case=%d e=%d got=%La want=%La\n", t, e, got, want);
            prv("input", v);
            return false;
        }
    }
    return true;
}

// 暴力对拍 + 退化/大规模用例;通过则 PASSED("最近点对")。
int main() {
    printf("== 最近点对.check:O(n^2) 暴力对拍 + 退化/大规模用例 ==\n");
    CHECK(exponent_ranges(), "独立hypot参考3000组：全指数/平方下溢边界/大平移/远近混合");

    // ===== 1. 解析用例 =====
    {
        CHECK(dis_closest({P(3, 4)}) > 1e99L, "dis_closest:n=1 -> 返回 1e100 的无解标记");
        CHECK(dis_closest({}) > 1e99L, "dis_closest:空输入 -> 无解标记");
        CHECK(eqd(dis_closest({P(0, 0), P(3, 4)}), 5, 1e-12L), "dis_closest:n=2 -> 两点距离");
        CHECK(eqd(dis_closest({P(0, 0), P(3, 4), P(0, 0)}), 0, 1e-12L), "dis_closest:有重复点 -> 0");
        CHECK(eqd(dis_closest({P(0, 0), P(1, 1), P(2, 2), P(10, 10)}), sqrtl(2.0L), 1e-12L), "dis_closest:全共线 -> 相邻两点");
        CHECK(eqd(dis_closest({P(0, 0), P(1000000000, 0), P(0, 1000000000), P(1000000000, 1000000000), P(999999999, 1000000000)}),
                  1, 1e-12L), "dis_closest:±1e9 端点 + 一个 1 距离的点 -> 1");
        {
            vector<p2> g;
            For(i, 0, 9) For(j, 0, 9) g.push_back(P(i * 7, j * 7));
            CHECK(eqd(dis_closest(g), 7, 1e-12L), "dis_closest:10×10 正方形网格(间距 7)-> 7");
            g.push_back(P(3.5L, 3.5L));
            CHECK(eqd(dis_closest(g), dis(P(3.5L, 3.5L) - P(0, 0)), 1e-12L), "dis_closest:网格 + 格中心点 -> 中心到角点的距离");
        }
    }

    // ===== 2. 小整数坐标(大量重复)与暴力对拍 =====
    {
        int bad = 0;
        For(t, 1, 4000) {
            int n = (int) rnd(2, 40), R = (int) rnd(1, 8);
            vector<p2> v;
            if(t % 4 == 0) {  // 每 4 组里有一组「保证互异」
                for(int i = 0; i < n; ++i) v.push_back(P((db) rnd(-200, 200), (db) rnd(-200, 200)));
                for(int i = 0; i < n; ++i) for(int j = 0; j < i; ++j) if(v[i] == v[j]) v[i].x += 1000;
            } else {
                for(int i = 0; i < n; ++i) v.push_back(P((db) rnd(-R, R), (db) rnd(-R, R)));
            }
            db got = dis_closest(v), want = bruteClosest(v);
            if(!eqd(got, want, 1e-9 * (1 + want))) {
                if(bad < 3) printf("      [反例] 第 %d 组:模板 %.17Lg,暴力 %.17Lg\n", t, got, want), prv("输入", v);
                ++bad;
            }
        }
        CHECK(bad == 0, "dis_closest 与 O(n²) 暴力一致(4000 组:n=2..40,坐标 ±8 大量重复点 + 每 4 组一组互异点)");
    }

    // ===== 3. 中等规模随机分布与暴力对拍 =====
    {
        int bad = 0;
        For(t, 1, 3000) {
            int n = (int) rnd(2, 150);
            vector<p2> v;
            int kind = (int) rnd(0, 2);
            if(kind == 0) {
                for(int i = 0; i < n; ++i) v.push_back(P((db) rnd(-1000000, 1000000), (db) rnd(-1000000, 1000000)));
            } else if(kind == 1) {  // 圆上点:最近点对常常落在「相邻弧」上
                db R = (db) rnd(1, 10000);
                for(int i = 0; i < n; ++i) v.push_back(P(R * cos(2 * pi * i / n), R * sin(2 * pi * i / n)));
            } else {  // 细长分布:横坐标差很小、纵坐标差很大
                for(int i = 0; i < n; ++i) v.push_back(P((db) rnd(-5, 5), (db) rnd(-1000000, 1000000)));
            }
            db got = dis_closest(v), want = bruteClosest(v);
            if(!eqd(got, want, 1e-9 * (1 + want))) {
                if(bad < 3) printf("      [反例] 第 %d 组(kind=%d,n=%d):模板 %.17Lg,暴力 %.17Lg\n", t, kind, n, got, want), prv("输入", v);
                ++bad;
            }
        }
        CHECK(bad == 0, "dis_closest 与 O(n²) 暴力一致(3000 组:n=2..150,随机/圆上点/细长分布三类,相对 1e-9)");
    }

    // ===== 4. 性质级断言 =====
    {
        int bad = 0;
        For(t, 1, 3000) {
            int n = (int) rnd(2, 60);
            vector<p2> v;
            for(int i = 0; i < n; ++i) v.push_back(P((db) rnd(-1000, 1000), (db) rnd(-1000, 1000)));
            db got = dis_closest(v);
            if(got < 0) ++bad;
            for(int i = 0; i < min(n, 25); ++i) for(int j = i + 1; j < min(n, 25); ++j) {  // 答案必须 <= 任意一对具体点的距离
                db d = dis(v[i] - v[j]);
                if(got > d + 1e-9 * (1 + d)) {
                    if(bad < 3) printf("      [反例] 答案 %.17Lg 大于点对距离 %.17Lg\n", got, d), prv("输入", v);
                    ++bad;
                }
            }
            bool dup = false;  // 答案 == 0 <=> 存在重复点
            for(int i = 0; i < n; ++i) for(int j = 0; j < i; ++j) if(v[i] == v[j]) dup = true;
            if(dup != (got == 0)) {
                if(bad < 3) printf("      [反例] 有重复点=%d,而答案=%Lg\n", (int) dup, got), prv("输入", v);
                ++bad;
            }
        }
        CHECK(bad == 0, "dis_closest 性质:非负、<= 任意点对距离、== 0 当且仅当有重复点(3000 组)");
    }

    // ===== 5. 大规模:真答案已知的晶格 + 一对跨中位线的近距离点 =====
    {
        vector<p2> v;
        int nx = 450, ny = 450;  // 202500 个晶格点
        for(int i = 0; i < nx; ++i) for(int j = 0; j < ny; ++j) v.push_back(P((db) i * 10000, (db) j * 10000));
        db cx = (db) (nx - 1) * 10000 / 2, cy = (db) (ny - 1) * 10000 / 2 + 5000;  // 落在两行晶格正中
        v.push_back(P(cx - 1.75L, cy)), v.push_back(P(cx + 1.75L, cy));
        db got = dis_closest(v);
        CHECK(eqd(got, 3.5L, 1e-9L), "dis_closest:20 万晶格点 + 一对跨中位线、距离 3.5 的点 -> 必须恰好是 3.5(大规模 + 跨分界线窄条)");
        vector<p2> samep(200000, P(-1234.5L, 9876.25L));
        CHECK(eqd(dis_closest(samep), 0, 1e-12L), "dis_closest:20 万个完全相同的点 -> 0(不许退化成 O(n²))");
        vector<p2> line;  // 20 万个共线点(横坐标递增很慢),最近对是相邻两点
        for(int i = 0; i < 200000; ++i) line.push_back(P((db) i * 0.001L, 1));
        CHECK(eqd(dis_closest(line), 0.001L, 1e-9L), "dis_closest:20 万个几乎共线的点 -> 相邻两点距离 0.001");
    }

    // ===== 6. ±1e9 大坐标 =====
    {
        int bad = 0;
        For(t, 1, 2000) {
            int n = (int) rnd(2, 60);
            vector<p2> v;
            for(int i = 0; i < n; ++i) v.push_back(P((db) rnd(-1000000000, 1000000000), (db) rnd(-1000000000, 1000000000)));
            if(t % 5 == 0) v.push_back(v[0]);  // 掺重复点
            db got = dis_closest(v), want = bruteClosest(v);
            if(!eqd(got, want, 1e-9 * (1 + want))) {
                if(bad < 3) printf("      [反例] 大坐标第 %d 组:模板 %.17Lg,暴力 %.17Lg\n", t, got, want), prv("输入", v);
                ++bad;
            }
        }
        CHECK(bad == 0, "dis_closest:±1e9 大坐标与暴力一致(2000 组:随机点 + 掺重复点,相对 1e-9)");
    }
    PASSED("最近点对");
}
