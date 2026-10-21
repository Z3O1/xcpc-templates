// init()/solveg()/f()/solve(): 数论 · Min_25 筛

// 求积性函数在 1..n 的前缀和，1 <= n <= 1e10；solve 可重复调用。
constexpr int i6 = (M + 1) / 6;
void sub(int &a, int b) { a -= b, a < 0 && (a += M); }
const int N = 1e5 + 10;
int p[N], pc, p2[N], s1[N], s2[N], s[N], g1[N * 2], g2[N * 2];
ll lim, n0;
void init() {
    if(pc) return;
    static const int n = 1e5;
    static int b[N];
    For(i, 2, n) {
        if(!b[i]) p[++pc] = i;
        for(int j = 1, x; j <= pc && (x = i * p[j]) <= n; ++j) {
            b[x] = 1;
            if(i % p[j] == 0) break;
        }
    }
    p[pc + 1] = n + 1;
    For(i, 1, pc) {
        p2[i] = 1ll * p[i] * p[i] % M;
        s1[i] = (s1[i - 1] + p[i]) % M, s2[i] = (s2[i - 1] + p2[i]) % M;
        s[i] = (s2[i] - s1[i] + M) % M;
    }
}
int id(ll n) { return n <= lim ? n : n0 / n + lim; }
void solveg(ll _n0) {
    n0 = _n0, lim = sqrt(n0);
    bool fl = n0 / lim == lim;
    auto put = [&](ll k, int slot) {
        k %= M;
        sub(g1[slot] = k * (k + 1) / 2 % M, 1);
        sub(g2[slot] = k * (k + 1) % M * (k + k + 1) % M * i6 % M, 1);
    };
    auto doit = [&](int m, int i1, int i2) {
        sub(g1[i1], 1ll * (g1[i2] - s1[m - 1] + M) * p[m] % M);
        sub(g2[i1], 1ll * (g2[i2] - s2[m - 1] + M) * p2[m] % M);
    };
    For(i, 1, lim) put(i, i);
    For(i, 1, lim - fl) put(n0 / i, i + lim);
    For(i, 1, pc) {
        ll x = p[i], y = x * x, z = n0 / x;
        if(y > n0) break;
        For(j, 1, min(lim - fl, n0 / y)) doit(i, j + lim, id(z / j));
        if(y <= lim) rFor(j, lim, y) doit(i, j, j / x);
    }
}
int g(ll n) { return (g2[id(n)] - g1[id(n)] + M) % M; }
ll func(ll n) { return n %= M, n * (n - 1) % M; }
int f(ll n, int m) {
    if(p[m] >= n) return 0;
    ll ans = g(n) - s[m] + M;
    For(j, m + 1, pc) {
        int q = p[j];
        if(1ll * q * q > n) break;
        for(ll w = q; w * q <= n; w *= q) ans += func(w) * f(n / w, j) % M + func(w * q);
    }
    return ans % M;
}
ll solve(ll n) { return init(), solveg(n), (f(n, 0) + 1) % M; }
