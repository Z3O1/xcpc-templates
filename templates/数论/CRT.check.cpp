// CRT 自测:与暴力枚举对照(含无解、非互质模、极值)
#include "../_check_base.hpp"
#include "exgcd.cpp"
#include "CRT.cpp"

int main() {
    // 小范围穷举:有解时必须等于 [0, lcm) 内最小解;无解必须都判无解
    int ok_cnt = 0, none_cnt = 0;
    For(p1, 1, 40) For(p2, 1, 40) For(a1, 0, p1 - 1) For(a2, 0, p2 - 1) {
        ll r = CRT(a1, p1, a2, p2);
        ll lcm = (ll) p1 / std::gcd(p1, p2) * p2, want = -1;
        For(x, 0, lcm - 1) if(x % p1 == a1 && x % p2 == a2) { want = x; break; }
        if(want < 0) {
            if(r != -1) { printf("  [FAIL] 应无解 a1=%d p1=%d a2=%d p2=%d got=%lld\n", a1, p1, a2, p2, r); return 1; }
            ++none_cnt;
        } else {
            if(r != want) { printf("  [FAIL] a1=%d p1=%d a2=%d p2=%d want=%lld got=%lld\n", a1, p1, a2, p2, want, r); return 1; }
            ++ok_cnt;
        }
    }
    printf("  [ok] 穷举 p1,p2<=40:有解 %d 组、无解 %d 组\n", ok_cnt, none_cnt);
    if(!ok_cnt || !none_cnt) return printf("  [FAIL] 用例覆盖不足\n"), 1;

    // 结果必须落在 [0, lcm)
    For(t, 1, 50000) {
        ll p1 = rnd(1, (ll) 1e9), p2 = rnd(1, (ll) 1e9);
        ll a1 = rnd(0, p1 - 1), a2 = rnd(0, p2 - 1);
        ll r = CRT(a1, p1, a2, p2);
        if(r < 0) continue;
        ll lcm = p1 / std::gcd(p1, p2) * p2;
        if(r < 0 || r >= lcm || r % p1 != a1 % p1 || r % p2 != a2 % p2)
            return printf("  [FAIL] 随机 p1=%lld p2=%lld r=%lld lcm=%lld\n", p1, p2, r, lcm), 1;
    }
    ok("5 万组随机大模数:落在 [0,lcm) 且同余成立");

    // 两两互质(经典 CRT)与极端边界
    {
        ll r = CRT(2, 3, 3, 5);
        CHECK(r == 8, "CRT(2 mod 3, 3 mod 5) = 8");
        r = CRT(1, 2, 2, 4);
        CHECK(r == -1, "CRT(x≡1 mod2, x≡2 mod4) 无解");
        r = CRT(0, 1, 0, 1);
        CHECK(r == 0, "模 1 退化");
        r = CRT(4, 6, 10, 15);
        CHECK(r == 10, "非互质有解 CRT(x≡4 mod6, x≡10 mod15) = 10");
    }

    // 多方程串联(文档里的用法),与暴力对照
    For(t, 1, 20000) {
        int n = (int) rnd(2, 4);
        vect<ll> a(n), p(n);
        ll mod = 1;
        For(i, 0, n - 1) p[i] = rnd(1, 20), a[i] = rnd(0, p[i] - 1), mod = mod / std::gcd(mod, p[i]) * p[i];
        if(mod > 20000) continue;
        ll ans = 0, m = 1;
        bool bad = 0;
        For(i, 0, n - 1) {
            ans = CRT(ans, m, a[i], p[i]);
            if(ans < 0) { bad = 1; break; }
            m = m / std::gcd(m, p[i]) * p[i];
        }
        ll want = -1;
        For(x, 0, m - 1) {
            bool good = 1;
            For(i, 0, n - 1) good &= x % p[i] == a[i];
            if(good) { want = x; break; }
        }
        if(bad) { if(want >= 0) return printf("  [FAIL] 误判无解\n"), 1; }
        else if(ans != want) return printf("  [FAIL] 多方程 ans=%lld want=%lld\n", ans, want), 1;
    }
    ok("多方程串联(2~4 个方程)与暴力一致");

    PASSED("CRT");
}
