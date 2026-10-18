// 多项式复合逆自测:双向独立复合、逐项反解/二次闭式、随机/长度/缩放、共享 NTT 集成
#include "../_check_base.hpp"
#include "ntt.cpp"
#include "多项式复合.cpp"
#include "多项式复合逆.cpp"

static void dump(const char *name, const poly &a) {
    printf("%s[%zu] =", name, a.size());
    for(int x : a) printf(" %d", x);
    puts("");
}
static void require_same(const char *what, const poly &f, int n, const poly &want, const poly &got) {
    if(want == got) return;
    printf("[FAIL] %s n=%d\n", what, n);
    dump("f", f);
    dump("want", want);
    dump("got", got);
    exit(1);
}
static poly ser_x(int n) {
    poly x(n);
    if(n > 1) x[1] = 1;
    return x;
}
// Horner 独立参考；不使用 NTT 或 mul。
static poly naive_comp(const poly &f, const poly &g, int n) {
    poly r(n);
    int last = min(n, (int)f.size()) - 1;
    while(last >= 0 && !f[last]) --last;
    rFor(i, last, 0) {
        poly t(n);
        For(a, 0, n - 1) if(r[a])
            For(b, 0, min(n - a, (int)g.size()) - 1) if(g[b])
                t[a + b] = (t[a + b] + (ll)r[a] * g[b]) % M;
        t[0] = (t[0] + f[i]) % M;
        r.swap(t);
    }
    return r;
}
// 第二种独立复合参考：累加 F[i]G^i。
static poly naive_comp2(const poly &f, const poly &g, int n) {
    poly r(n), p(n);
    p[0] = 1;
    For(i, 0, min(n, (int)f.size()) - 1) {
        For(j, 0, n - 1) r[j] = (r[j] + (ll)f[i] * p[j]) % M;
        poly t(n);
        For(a, 0, n - 1) if(p[a])
            For(b, 0, min(n - a, (int)g.size()) - 1) if(g[b])
                t[a + b] = (t[a + b] + (ll)p[a] * g[b]) % M;
        p.swap(t);
    }
    return r;
}
// 独立逐项反解：第 k 项仅以 f[1]*g[k] 线性出现。
static poly naive_inverse(const poly &f, int n) {
    poly g(n);
    int inv = qpow(f[1], M - 2, M);
    if(n > 1) g[1] = inv;
    For(k, 2, n - 1) {
        int c = naive_comp(f, g, k + 1)[k];
        g[k] = (ll)(M - c) * inv % M;
    }
    return g;
}
static poly rnd_f(int n, int mode) {
    poly f(max(n, 2));
    For(i, 1, n - 1) {
        f[i] = rnd(0, M - 1);
        if(mode == 1 && (i & 1)) f[i] = 0;
        if(mode == 2) f[i] = rnd(0, 4);
        if(mode == 3) f[i] = M - 1;
        if(mode == 4) f[i] = 0;
    }
    f[1] = mode == 2 ? 1 : rnd(1, M - 1);
    return f;
}
static void compare(const poly &f, int n, bool reverse_way) {
    poly saved = f;
    poly g = comp_inv(f, n);
    require_same("不改输入", f, n, saved, f);
    require_same("f(g)=x", f, n, ser_x(n), naive_comp(f, g, n));
    if(reverse_way) require_same("g(f)=x", f, n, ser_x(n), naive_comp(g, f, n));
    if((int)g.size() != n || g[0] != 0 || (n > 1 && g[1] != qpow(f[1], M - 2, M))) {
        dump("f", f);
        dump("g", g);
        exit(1);
    }
    for(int x : g)
        if(x < 0 || x >= M) {
            dump("g", g);
            exit(1);
        }
}

int main() {
    rng.seed(2026040102);
    require_same("参考手算", {0, 1, 1}, 4, {0, 1, 0, 0}, naive_comp({0, 1, 1}, {0, 1, M - 1, 2}, 4));
    For(t, 1, 2000) {
        int n = rnd(2, 14);
        poly f = rnd_f(n, t % 4), g = rnd_f(n, (t + 1) % 4);
        require_same("两个复合参考互校", f, n, naive_comp2(f, g, n), naive_comp(f, g, n));
        if(t <= 400) require_same("逐项反解参考", f, n, naive_inverse(f, n), comp_inv(f, n));
    }
    For(t, 1, 4000) {
        int n = rnd(2, 30);
        poly f = rnd_f(n, t % 4);
        if(t % 3 == 0) f[1] = 1;
        compare(f, n, true);
    }
    For(t, 1, 600) {
        int n = t <= 400 ? rnd(31, 60) : rnd(61, 100);
        poly f = rnd_f(n, t % 4);
        if(t % 2 == 0) f[1] = 1;
        compare(f, n, t <= 400 && t % 5 == 0);
    }
    For(t, 1, 22) {
        int n = t <= 20 ? rnd(150, 300) : rnd(400, 512);
        compare(rnd_f(n, t % 2), n, false);
    }
    ok("双参考互校、400 组逐项反解及 4622 组随机双向/正向独立复合验证");
    // 幂投影单独对拍：逐次朴素相乘生成 f^k，直接提取目标系数。
    For(t, 1, 1000) {
        int n = rnd(1, 35);
        poly f = rnd_f(n, t % 5);
        f.resize(n);
        poly power(n), want(n);
        power[0] = 1;
        For(k, 0, n - 1) {
            want[k] = power[n - 1];
            poly next(n);
            For(i, 0, n - 1)
                For(j, 0, n - 1 - i) next[i + j] = (next[i + j] + (ll)power[i] * f[j]) % M;
            power.swap(next);
        }
        require_same("幂投影逐系数独立参考", f, n, want, reversion::pw_pj(f, n));
    }
    ok("1000 组幂投影独立朴素参考，含 n=1、奇偶长度和稀疏输入");

    for(int n : {1, 2, 3, 6, 8, 17, 32, 33}) {
        for(int c : {1, 3, M - 1}) {
            poly want(n);
            if(n > 1) want[1] = qpow(c, M - 2, M);
            require_same("c*x 的闭式", {0, c}, n, want, comp_inv({0, c}, n));
            compare({0, c}, n, true);
        }
        poly f(n, 1), want(n);
        f[0] = 0;
        if(n == 1) f.push_back(1);
        For(i, 1, n - 1) want[i] = i & 1 ? 1 : M - 1;
        require_same("x/(1-x) 的逆", f, n, want, comp_inv(f, n));
    }
    require_same("x+x² 手算", {0, 1, 1}, 4, {0, 1, M - 1, 2}, comp_inv({0, 1, 1}, 4));
    require_same("负几何级数自反", {0, M - 1, M - 1, M - 1, M - 1}, 5,
                 {0, M - 1, M - 1, M - 1, M - 1}, comp_inv({0, M - 1, M - 1, M - 1, M - 1}, 5));
    // 保留旧测试的契约外零一次项观察，不把它解释成存在复合逆。
    for(int n : {2, 3, 8, 50, 200}) {
        poly f = rnd_f(n, 0);
        f[1] = 0;
        require_same("契约外 F[1]=0 的既有行为", f, n, poly(n), comp_inv(f, n));
    }
    ok("n=1/2、短输入、线性/几何级数/二次手算和既有契约外观察");

    for(int n : {500, 1000, 2000}) {
        poly f(8);
        f[1] = 1, f[3] = rnd(1, M - 1), f[7] = rnd(1, M - 1);
        compare(f, n, false);
    }
    {
        int n = 2000;
        poly want(n);
        want[1] = 1;
        For(k, 2, n - 1) {
            ll s = 0;
            For(i, 1, k - 1) s += (ll)want[i] * want[k - i] % M;
            want[k] = (M - s % M) % M;
        }
        require_same("g+g²=x 的独立二次递推", {0, 1, 1}, n, want, comp_inv({0, 1, 1}, n));
    }
    for(int n : {1024, 2048}) {
        poly f = rnd_f(n, 0), scaled(n);
        int inv_c = qpow(13, M - 2, M);
        ll pw = 1;
        For(i, 1, n - 1) scaled[i] = (ll)f[i] * pw % M, pw = pw * inv_c % M;
        poly want = comp_inv(f, n);
        pw = 13;
        For(i, 0, n - 1) want[i] = want[i] * pw % M, pw = pw * inv_c % M;
        require_same("c*F(x/c) 的缩放对称性", scaled, n, want, comp_inv(scaled, n));
    }
    For(t, 1, 200) {
        int n = rnd(2, 40);
        poly f = rnd_f(n, t % 4), long_f = f;
        long_f.resize(n + 5);
        For(i, n, n + 4) long_f[i] = rnd(0, M - 1);
        poly want = comp_inv(f, n);
        require_same("重复调用", f, n, want, comp_inv(f, n));
        require_same("高次不倒入低次", long_f, n, want, comp_inv(long_f, n));
        f.resize(rnd(2, n));
        poly padded = f;
        padded.resize(n);
        require_same("短输入补零", f, n, comp_inv(padded, n), comp_inv(f, n));
    }
    for(const poly &f : {poly{0, 1, 0}, poly{0, 1, 1, 0, 0}, poly{0, 1, 0, 0, 7}}) {
        int n = f.size() - 1;
        poly truncated = f;
        truncated.resize(n);
        require_same("长输入历史反例", f, n, comp_inv(truncated, n), comp_inv(f, n));
    }
    ok("大 n 稀疏独立复合、二次递推、缩放对称性及长短输入硬回归");

    // 同一单元同时使用两份 FPS 模板和普通乘法；覆盖两侧分治边界及先大后小。
    for(int n : {8193, 8192, 8191, 1025, 1024, 1023, 33, 2, 1, 2049}) {
        poly f = rnd_f(n, 0), g = comp_inv(f, n);
        require_same("共享 NTT: f(g)=x", f, n, ser_x(n), comp(f, g, n));
        require_same("共享 NTT: g(f)=x", f, n, ser_x(n), comp(g, f, n));
    }
    {
        poly warm = mul(poly(5, 1), poly(5, 1));
        require_same("大乘法后的常数乘法", {1}, 1, {1}, mul({1}, {1}));
        For(t, 1, 2400) {
            poly a(rnd(0, 70)), b(rnd(0, 70));
            for(int &x : a) x = rnd(0, M - 1);
            for(int &x : b) x = rnd(0, M - 1);
            poly want;
            if(!a.empty() && !b.empty()) {
                want.resize(a.size() + b.size() - 1);
                For(i, 0, (int)a.size() - 1)
                    For(j, 0, (int)b.size() - 1) want[i + j] = (want[i + j] + (ll)a[i] * b[j]) % M;
            }
            if(want != mul(a, b)) dump("b", b);
            require_same("公共乘法独立参考", a, want.size(), want, mul(a, b));
        }
    }
    ok("双模板集成、先大后小共享缓存与 2400 组公共乘法独立参考");
    // F=x/(1-x) 的逆为 x/(1+x)，旧版上限处仍逐项核对而非只比较哈希。
    for(int n : {524288, 131073, 3}) {
        poly f(n, 1), want(n);
        f[0] = 0;
        For(i, 1, n - 1) want[i] = i & 1 ? 1 : M - 1;
        require_same("旧版容量的稠密几何级数复合逆", f, n, want, comp_inv(f, n));
    }
    ok("2^22 内部变换容量、奇数大长度、稠密输入独立闭式");
    PASSED("多项式复合逆");
}
