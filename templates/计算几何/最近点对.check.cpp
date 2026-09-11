// 最近点对 自测:dis_closest 与 O(n²) 全点对暴力对拍 + 退化/极端/大规模用例
//
// 对拍方式(文件内自写,不用模板):
//   bruteClosest(v):O(n²) 双重循环逐个算 dis() 取最小(独立写一遍,连排序都不用)。
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
#include "../_check_base.hpp"
#include "geo.cpp"
#include "最近点对.cpp"

static p2 P(db x, db y) { return {x, y}; }
static db Abs(db x) { return x < 0 ? -x : x; }
static bool eqd(db a, db b, db tol) { return Abs(a - b) < tol; }
static void pr(p2 a) { printf("(%g,%g)", (double) a.x, (double) a.y); }
static void prv(const char *s, const vector<p2> &v) {
    printf("      %s(%zu):", s, v.size());
    ForD(i, 0, (int) min<size_t>(v.size(), 20)) printf(" "), pr(v[i]);
    printf("%s\n", v.size() > 20 ? " ..." : "");
}

// ---------- 独立参考实现:O(n²) 暴力 ----------
static db bruteClosest(const vector<p2> &v) {
    int n = v.size();
    if(n < 2) return -1;
    db ans = 1e100L;
    ForD(i, 0, n) ForD(j, i + 1, n) ans = min(ans, dis(v[i] - v[j]));
    return ans;
}

int main() {
    printf("== 最近点对.check:O(n^2) 暴力对拍 + 退化/大规模用例 ==\n");

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
                ForD(i, 0, n) v.push_back(P((db) rnd(-200, 200), (db) rnd(-200, 200)));
                ForD(i, 0, n) ForD(j, 0, i) if(v[i] == v[j]) v[i].x += 1000;
            } else {
                ForD(i, 0, n) v.push_back(P((db) rnd(-R, R), (db) rnd(-R, R)));
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
                ForD(i, 0, n) v.push_back(P((db) rnd(-1000000, 1000000), (db) rnd(-1000000, 1000000)));
            } else if(kind == 1) {  // 圆上点:最近点对常常落在「相邻弧」上
                db R = (db) rnd(1, 10000);
                ForD(i, 0, n) v.push_back(P(R * cos(2 * pi * i / n), R * sin(2 * pi * i / n)));
            } else {  // 细长分布:横坐标差很小、纵坐标差很大
                ForD(i, 0, n) v.push_back(P((db) rnd(-5, 5), (db) rnd(-1000000, 1000000)));
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
            ForD(i, 0, n) v.push_back(P((db) rnd(-1000, 1000), (db) rnd(-1000, 1000)));
            db got = dis_closest(v);
            if(got < 0) ++bad;
            ForD(i, 0, min(n, 25)) ForD(j, i + 1, min(n, 25)) {  // 答案必须 <= 任意一对具体点的距离
                db d = dis(v[i] - v[j]);
                if(got > d + 1e-9 * (1 + d)) {
                    if(bad < 3) printf("      [反例] 答案 %.17Lg 大于点对距离 %.17Lg\n", got, d), prv("输入", v);
                    ++bad;
                }
            }
            bool dup = false;  // 答案 == 0 <=> 存在重复点
            ForD(i, 0, n) ForD(j, 0, i) if(v[i] == v[j]) dup = true;
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
        ForD(i, 0, nx) ForD(j, 0, ny) v.push_back(P((db) i * 10000, (db) j * 10000));
        db cx = (db) (nx - 1) * 10000 / 2, cy = (db) (ny - 1) * 10000 / 2 + 5000;  // 落在两行晶格正中
        v.push_back(P(cx - 1.75L, cy)), v.push_back(P(cx + 1.75L, cy));
        db got = dis_closest(v);
        CHECK(eqd(got, 3.5L, 1e-9L), "dis_closest:20 万晶格点 + 一对跨中位线、距离 3.5 的点 -> 必须恰好是 3.5(大规模 + 跨分界线窄条)");
        vector<p2> samep(200000, P(-1234.5L, 9876.25L));
        CHECK(eqd(dis_closest(samep), 0, 1e-12L), "dis_closest:20 万个完全相同的点 -> 0(不许退化成 O(n²))");
        vector<p2> line;  // 20 万个共线点(横坐标递增很慢),最近对是相邻两点
        ForD(i, 0, 200000) line.push_back(P((db) i * 0.001L, 1));
        CHECK(eqd(dis_closest(line), 0.001L, 1e-9L), "dis_closest:20 万个几乎共线的点 -> 相邻两点距离 0.001");
    }

    // ===== 6. ±1e9 大坐标 =====
    {
        int bad = 0;
        For(t, 1, 2000) {
            int n = (int) rnd(2, 60);
            vector<p2> v;
            ForD(i, 0, n) v.push_back(P((db) rnd(-1000000000, 1000000000), (db) rnd(-1000000000, 1000000000)));
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
