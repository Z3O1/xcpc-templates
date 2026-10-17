// FFT(复数版卷积)自测:与 O(nm) 暴力卷积对拍 + 精度边界
#include "../_check_base.hpp"
#include "FFT.cpp"

// 运行 FFT 卷积自测；覆盖随机对拍、长度/空输入/精度边界；失败返回非零。
int main() {
    // 1) 随机对拍(系数小,double 精度足够)
    For(t, 1, 30000) {
        int n = (int) rnd(1, 40), m = (int) rnd(1, 40);
        vector<db> a(n, 0), b(m, 0);
        For(i, 0, n - 1) a[i] = (db) rnd(-1000, 1000);
        For(i, 0, m - 1) b[i] = (db) rnd(-1000, 1000);
        vector<db> got = fft_conv(a, b), want(n + m - 1, 0);
        For(i, 0, n - 1) For(j, 0, m - 1) want[i + j] += a[i] * b[j];
        if((int) got.size() != n + m - 1) return printf("  [FAIL] 长度 %zu 应 %d\n", got.size(), n + m - 1), 1;
        For(i, 0, n + m - 2) {
            db err = fabs(got[i] - want[i]) / max((db) 1, fabs(want[i]));
            if(!isfinite(got[i]) || err > 1e-9) return printf("  [FAIL] t=%d n=%d m=%d 第 %d 项 %.10Lf 应 %.10Lf\n", t, n, m, i, got[i], want[i]), 1;
        }
    }
    ok("3 万组随机小系数卷积(相对误差 < 1e-9)");

    // 2) 长度边界:直接FFT空/长度1、1×1、1×k、k×1、非2的幂卷积长度
    {
        vector<cp> single{cp(3, -2)}, empty;
        fft(single, 1), fft(single, -1);
        CHECK(single.size() == 1 && single[0] == cp(3, -2), "长度1的正逆FFT均为恒等");
        fft(empty, 1), fft(empty, -1);
        CHECK(empty.empty(), "缓存曾使用大长度后,空FFT仍安全");
        vector<db> a{3}, b{4};
        CHECK(fft_conv(a, b).size() == 1 && fabs(fft_conv(a, b)[0] - 12) < 1e-9, "1×1 卷积 = 12");
        for(int m : {1, 2, 3, 5, 17, 33, 100}) {
            vector<db> x{2}, y(m, 1);
            vector<db> g = fft_conv(x, y);
            if((int) g.size() != m) return printf("  [FAIL] 1×%d 长度\n", m), 1;
            For(i, 0, m - 1) if(!isfinite(g[i]) || fabs(g[i] - 2) > 1e-9) return printf("  [FAIL] 1×%d 值\n", m), 1;
        }
        ok("1×1 / 1×k / 非 2 的幂长度");
    }

    // 3) 空输入与含 0 元素
    CHECK(fft_conv({}, {1}).empty() && fft_conv({1}, {}).empty(), "空输入返回空");
    {
        vector<db> a{0, 0, 0}, b{1, 2};
        vector<db> g = fft_conv(a, b);
        if((int) g.size() != 4) return printf("  [FAIL] 全零长度\n"), 1;
        For(i, 0, 3) if(!isfinite(g[i]) || fabs(g[i]) > 1e-9) return printf("  [FAIL] 全零结果非零\n"), 1;
        ok("空输入 + 全零多项式");
    }

    // 4) 精度边界:系数 ~1e6、长度 1024(结果 ~1e15,仍在 double 能力内但误差会放大)
    {
        int n = 1024;
        vector<db> a(n, 0), b(n, 0);
        For(i, 0, n - 1) a[i] = (db) rnd(-1000000, 1000000);
        For(i, 0, n - 1) b[i] = (db) rnd(-1000000, 1000000);
        vector<db> got = fft_conv(a, b);
        CHECK((int) got.size() == 2 * n - 1, "大卷积长度");
        for(db x : got) CHECK(isfinite(x), "有限输入的大卷积不得返回 NaN/Inf");
        // 抽样几个位置与暴力对照
        For(i, 0, 5) {
            int idx = (int) rnd(0, 2 * n - 2);
            i128 want = 0;
            for(int j = max(0, idx - n + 1); j <= min(idx, n - 1); ++j) want += (i128) (ll) a[j] * (ll) b[idx - j];
            db err = fabs(got[idx] - (db) want) / max((db) 1, fabs((db) want));
            if(!isfinite(got[idx]) || err > 1e-6) return printf("  [FAIL] 大系数 idx=%d got=%.3Lf want=%.3Lf\n", idx, got[idx], (db) want), 1;
        }
        ok("1024 长度、系数 1e6 量级(相对误差 < 1e-6)");
    }

    PASSED("FFT");
}
