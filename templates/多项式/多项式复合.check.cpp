// 多项式复合自测:双独立参考、随机/退化/代数性质、短输入与分治边界、共享 NTT 缓存
#include "../_check_base.hpp"
#include "ntt.cpp"
#include "多项式复合.cpp"

static void dump(const char *name, const poly &a) {
    printf("%s[%zu] =", name, a.size());
    for(int x : a) printf(" %d", x);
    puts("");
}
static void require_same(const char *what, const poly &f, const poly &g, int n, const poly &want, const poly &got) {
    if(want == got) return;
    printf("[FAIL] %s n=%d\n", what, n);
    dump("f", f);
    dump("g", g);
    dump("want", want);
    dump("got", got);
    exit(1);
}
// 独立 Horner 参考，不调用 mul；G[0] 按接口约定忽略。
static poly naive_comp(const poly &f, const poly &g, int n) {
    poly r(n);
    int last = min(n, (int)f.size()) - 1;
    while(last >= 0 && !f[last]) --last;
    rFor(i, last, 0) {
        poly t(n);
        For(a, 0, n - 1) if(r[a])
            For(b, 1, min(n - a, (int)g.size()) - 1) if(g[b])
                t[a + b] = (t[a + b] + (ll)r[a] * g[b]) % M;
        t[0] = (t[0] + f[i]) % M;
        r.swap(t);
    }
    return r;
}
// 第二种独立参考：逐次生成 G^i，再累加 F[i]G^i。
static poly naive_comp2(const poly &f, const poly &g, int n) {
    poly r(n), p(n);
    p[0] = 1;
    For(i, 0, min(n, (int)f.size()) - 1) {
        For(j, 0, n - 1) r[j] = (r[j] + (ll)f[i] * p[j]) % M;
        poly t(n);
        For(a, 0, n - 1) if(p[a])
            For(b, 1, min(n - a, (int)g.size()) - 1) if(g[b])
                t[a + b] = (t[a + b] + (ll)p[a] * g[b]) % M;
        p.swap(t);
    }
    return r;
}
static poly rnd_poly(int n, int mode) {
    poly a(n);
    For(i, 0, n - 1) {
        a[i] = rnd(0, M - 1);
        if(mode == 1 && (i & 1)) a[i] = 0;
        if(mode == 2) a[i] = 0;
        if(mode == 3) a[i] = M - 1;
        if(mode == 4) a[i] = rnd(0, 4);
    }
    if(mode == 5 && n) {
        fill(all(a), 0);
        a[rnd(0, n - 1)] = rnd(1, M - 1);
    }
    return a;
}
static poly rnd_g(int n, int mode) {
    poly g = rnd_poly(n, mode);
    if(n) g[0] = 0;
    return g;
}
static void compare(const poly &f, const poly &g, int n) {
    poly saved_f = f, saved_g = g;
    poly got = comp(f, g, n);
    require_same("Horner 对拍", f, g, n, naive_comp(f, g, n), got);
    require_same("不修改 f", f, g, n, saved_f, f);
    require_same("不修改 g", f, g, n, saved_g, g);
    for(int x : got)
        if(x < 0 || x >= M) {
            dump("got", got);
            exit(1);
        }
}

int main() {
    rng.seed(2026040101);
    require_same("参考手算", {1, 2, 3}, {0, 1, 1}, 5, {1, 2, 5, 6, 3}, naive_comp({1, 2, 3}, {0, 1, 1}, 5));
    For(t, 1, 2000) {
        int n = rnd(1, 16);
        poly f = rnd_poly(n, t % 6), g = rnd_g(n, (t + 3) % 6);
        require_same("两个参考互校", f, g, n, naive_comp2(f, g, n), naive_comp(f, g, n));
    }
    For(t, 1, 6000) {
        int n = rnd(1, 60);
        poly f = rnd_poly(n, t % 6), g = rnd_g(n, t / 6 % 6);
        if(t % 5 == 0) f.resize(rnd(0, n));
        if(t % 7 == 0) g.resize(rnd(0, n));
        compare(f, g, n);
    }
    For(t, 1, 600) {
        int n = rnd(61, 120);
        poly f = rnd_poly(n, t % 6), g = rnd_g(n, (t + 1) % 6);
        if(t % 4 == 0) f.resize(rnd(0, n));
        compare(f, g, n);
    }
    For(t, 1, 65) {
        int n = t <= 60 ? rnd(150, 300) : rnd(400, 512);
        compare(rnd_poly(n, t % 3), rnd_g(n, (t + 2) % 3), n);
    }
    ok("双参考互校及 6665 组随机逐系数对拍，含空/短 F、G");

    for(int n : {1, 2, 3, 4, 7, 8, 9, 31, 32, 33}) {
        for(poly f : {poly{}, poly{7}, poly{0, 1}, poly{0, 0, 1}, rnd_poly(n + 3, 3)}) {
            compare(f, {}, n);
            compare(f, {0}, n);
            compare(f, {0, 1}, n);
            compare(f, {0, 12345}, n);
            compare(f, {0, 0, 0, 1}, n);
        }
    }
    require_same("忽略 G[0]", {1, 2, 3}, {9, 1, 1}, 5, comp({1, 2, 3}, {0, 1, 1}, 5), comp({1, 2, 3}, {9, 1, 1}, 5));
    ok("退化、单项式、全 M-1、n=1 及非零 G[0] 的既有行为");

    For(t, 1, 200) {
        int n = rnd(2, 40);
        poly f = rnd_poly(n, t % 4), g = rnd_g(n, (t + 1) % 4), h = rnd_g(n, (t + 2) % 4);
        poly f2 = rnd_poly(n, (t + 3) % 4), sum(n);
        For(i, 0, n - 1) sum[i] = (f[i] + f2[i]) % M;
        poly a = comp(f, g, n), b = comp(f2, g, n);
        For(i, 0, n - 1) a[i] = (a[i] + b[i]) % M;
        require_same("对 F 线性", f, g, n, a, comp(sum, g, n));
        poly long_f = f, long_g = g;
        long_f.resize(2 * n);
        long_g.resize(2 * n);
        For(i, n, 2 * n - 1) long_f[i] = rnd(0, M - 1), long_g[i] = rnd(0, M - 1);
        require_same("高次忽略", long_f, long_g, n, comp(f, g, n), comp(long_f, long_g, n));
        poly big = comp(long_f, long_g, 2 * n);
        big.resize(n);
        require_same("截断一致", f, g, n, comp(f, g, n), big);
        require_same("结合律", f, g, n, comp(comp(f, g, n), h, n), comp(f, comp(g, h, n), n));
        require_same("右单位元", f, {0, 1}, n, f, comp(f, {0, 1}, n));
        require_same("左单位元", {0, 1}, g, n, g, comp({0, 1}, g, n));
        int k = rnd(0, min(n - 1, 6));
        poly monomial(k + 1);
        monomial[k] = 1;
        require_same("幂律", monomial, g, n, naive_comp2(monomial, g, n), comp(monomial, g, n));
        int c = 12345;
        poly want(n);
        ll pw = 1;
        For(i, 0, n - 1) want[i] = (ll)f[i] * pw % M, pw = pw * c % M;
        require_same("线性 G 闭式", f, {0, c}, n, want, comp(f, {0, c}, n));
    }
    ok("200 组线性、截断、高次忽略、结合律、单位元、幂律和缩放闭式");

    for(int n : {500, 1000, 2000}) {
        poly f(10), g = rnd_g(n, 0);
        f[1] = rnd(1, M - 1), f[3] = rnd(1, M - 1), f[9] = rnd(1, M - 1);
        compare(f, g, n);
    }
    // 按分治层数变化的两侧及先大后小顺序，校验稠密输入的 O(n) 闭式。
    for(int n : {65536, 4097, 4096, 4095, 1025, 1024, 1023, 17, 2, 1, 2049}) {
        poly f = rnd_poly(n + 3, 0), want(n);
        ll pw = 1;
        For(i, 0, n - 1) want[i] = (ll)f[i] * pw % M, pw = pw * 13 % M;
        require_same("大规模闭式及缓存变长", f, {0, 13}, n, want, comp(f, {0, 13}, n));
    }
    ok("大 n 独立参考、2 的幂两侧、先大后小共享 NTT 缓存");
    // F=1/(1-x), G=x/(1-x)：F(G)=(1-x)/(1-2x)，稠密双输入的独立闭式。
    // n=2^19 时内部 NTT 达到旧版 prep(22) 的 2^22 容量。
    for(int n : {524288, 131073, 3}) {
        poly f(n, 1), g(n, 1), want(n, 1);
        g[0] = 0;
        For(i, 2, n - 1) want[i] = 2ll * want[i - 1] % M;
        require_same("旧版容量的稠密几何级数复合", f, g, n, want, comp(f, g, n));
    }
    ok("2^22 内部变换容量、奇数大长度、稠密双输入独立闭式");
    PASSED("多项式复合");
}
