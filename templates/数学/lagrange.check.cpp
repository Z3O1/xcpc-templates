// lag()/lag_i()/lag_p() 自测:已知系数逐位还原、独立扩欧拉格朗日公式、整数点及负数归一边界
#include "../_check_base.hpp"

static const int FN = 4096;
static int FA[FN + 1], IFA[FN + 1];
// 查询预处理阶乘；要求 0 <= k <= FN；返回 k! mod M。
static int fac(int k) { return FA[k]; }
// 查询预处理逆阶乘；要求 0 <= k <= FN；返回 (k!)^-1 mod M。
static int ifac(int k) { return IFA[k]; }
// 预处理阶乘与逆阶乘；要求 M 为质数且 FN < M；覆盖 FA/IFA 全表。
static void init_fac() {
    FA[0] = 1;
    For(i, 1, FN) FA[i] = (ll)FA[i - 1] * i % M;
    IFA[FN] = qpow(FA[FN], M - 2);
    rFor(i, FN, 1) IFA[i - 1] = (ll)IFA[i] * i % M;
}
#include "lagrange.cpp"

// 代入求值；要求系数已归一到 [0,M),x 在模 M 意义下取值；返回多项式值。
static int evalp(const vector<int> &c, int x) {
    ll s = 0;
    rFor(i, (int)c.size() - 1, 0) s = (s * x + c[i]) % M;
    return s;
}
// 用扩欧求逆元；要求 v 与 M 互质；返回 [0,M) 内逆元。
static int ref_inv(int v) {
    ll a = v, b = M, x = 1, y = 0;
    while(b) {
        ll q = a / b;
        swap(a -= q * b, b), swap(x -= q * y, y);
    }
    return (x % M + M) % M;
}
// 独立拉格朗日求值；要求 a 为 0..n-1 点值且 n < M；返回 x 处值。
static int ref_lag_eval(const vector<int> &a, int x) {
    ll ans = 0;
    For(i, 0, (int)a.size() - 1) {
        ll num = a[i], den = 1;
        For(j, 0, (int)a.size() - 1) if(i != j) {
            num = num * (x - j + M) % M;
            den = den * (i - j + M) % M;
        }
        ans = (ans + num * ref_inv(den)) % M;
    }
    return ans;
}
// 生成随机多项式；要求 n >= 0；返回 n 项 [0,M) 系数。
static vector<int> random_poly(int n) {
    vector<int> c(n);
    for(int &x : c) x = rnd(0, M - 1);
    return c;
}
// 输出失败信息并终止；不返回,会打印相关多项式与期望/实际值。
static void fail(const char *what, const vector<int> &c, int x, int want, int got) {
    cerr << "[FAIL] " << what << " x/index=" << x << " want=" << want << " got=" << got << "\nc=";
    for(int v : c) cerr << ' ' << v;
    cerr << endl;
    exit(1);
}
// 校验 lag 还原结果；要求 pts 来自 c 且横坐标互异；失败时终止进程。
static void check_lag(const vector<int> &c, const vector<pair<int, int>> &pts) {
    vector<int> got = lag(pts);
    if(got.size() != pts.size() + 1) fail("lag 返回长度", c, 0, pts.size() + 1, got.size());
    For(i, 0, (int)c.size() - 1) if(got[i] != c[i]) fail("lag 系数", c, i, c[i], got[i]);
    if(got.back()) fail("lag 最高项", c, 0, 0, got.back());
    for(auto [x, y] : pts) if(evalp(got, x) != y) fail("lag 点值", c, x, y, evalp(got, x));
}
// 运行拉格朗日插值自测；覆盖任意点、连续点、原地插值、单点与空点约定。
int main() {
    rng.seed(20040924);
    init_fac();
    CHECK(fac(10) == 3628800 && (ll)fac(10) * ifac(10) % M == 1, "fac/ifac 自洽");
    For(t, 1, 20000) {
        int k = rnd(1, 9);
        vector<int> c = random_poly(k);
        vector<pair<int, int>> pts;
        while((int)pts.size() < k) {
            int x = (t % 3 == 0 && pts.empty()) ? 0 : rnd(0, M - 1);
            bool dup = false;
            for(auto p : pts) dup |= p.first == x;
            if(!dup) pts.push_back({x, evalp(c, x)});
        }
        check_lag(c, pts);
    }
    For(t, 1, 3000) {
        vector<int> c = random_poly(rnd(2, 20));
        vector<pair<int, int>> pts;
        For(i, 0, (int)c.size() - 1) pts.push_back({i, evalp(c, i)});
        check_lag(c, pts);
    }
    ok("2 万组任意互异点(含0) + 3000 组连续整数点,逐系数还原");
    For(t, 1, 20000) {
        int n = rnd(1, 40);
        vector<int> c = random_poly(n), a(n + 1);
        For(i, 1, n) a[i] = evalp(c, i);
        lag_i(n, a.data());
        For(i, 0, n - 1) if(a[i] != c[i]) fail("lag_i", c, i, c[i], a[i]);
        if(a[n]) fail("lag_i 末项清零", c, n, 0, a[n]);
    }
    ok("2 万组 lag_i:1..n 点值原地还原系数并清零末项");
    For(t, 1, 40000) {
        int n = rnd(1, 30), mode = t % 3;
        vector<int> c = random_poly(n), a(n);
        For(i, 0, n - 1) a[i] = evalp(c, i);
        int x = mode == 0 ? rnd(0, n - 1) : rnd(0, M - 1);
        if(mode == 2) x = (M - x) % M;
        int got = lag_p(n, a.data(), x), want = evalp(c, x);
        if(got != want) fail("lag_p", c, x, want, got);
        if(mode && ref_lag_eval(a, x) != got) fail("lag_p 独立公式", a, x, ref_lag_eval(a, x), got);
    }
    ok("4 万组 lag_p:命中整数点/任意点/负数归一,并与扩欧参考公式一致");
    CHECK(lag({{5, 9}}) == vector<int>({9, 0}), "单点(5,9)");
    CHECK(lag({{0, 7}}) == vector<int>({7, 0}), "单点横坐标0");
    CHECK(lag({}) == vector<int>({0}), "空点集的原有返回约定");
    vector<int> c{M - 1, M - 1, M - 1};
    check_lag(c, {{0, evalp(c, 0)}, {M - 1, evalp(c, M - 1)}, {M - 2, evalp(c, M - 2)}});
    vector<int> a{0, 9};
    lag_i(1, a.data());
    CHECK(a == vector<int>({9, 0}), "lag_i n=1");
    CHECK(lag_p(1, a.data(), 12345) == 9 && lag_p(1, a.data(), 0) == 9, "lag_p n=1");
    PASSED("lagrange");
}
