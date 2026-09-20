// init()/solveg()/f()/solve(): 数论 · Min_25 筛

const int M = 1e9 + 7;
constexpr int i2 = M + 1 >> 2, i6 = (M + 1) / 6;
void add(int &a, int b) {
    a += b;
    a >= M && (a -= M);
}
void sub(int &a, int b) {
    a -= b;
    a < 0 && (a += M);
}
const int N = 1e5 + 10;
int p[N], pc, p2[N];
int s1[N], s2[N], s[N];
void init() {
    static const int n = 1e5;
    static int b[N];
    For(i, 2, n) {
        if(!b[i]) p[++pc] = i;
        For(j, 1, pc) if(i * p[j] <= n) {
            b[i * p[j]] = 1;
            if(i % p[j] == 0) break;
        }
        else break;
    }
    p[pc + 1] = n + 1;
    For(i, 1, pc) {
        p2[i] = 1ll * p[i] * p[i] % M;
        s1[i] = (s1[i - 1] + p[i]) % M;
        s2[i] = (s2[i - 1] + p2[i]) % M;
        s[i] = (s2[i] - s1[i] + M) % M;
    }
}
ll lim, n0;
int g1[N * 2], g2[N * 2];
int id(ll n) {
    return n <= lim ? n : n0 / n + lim;
}
void solveg(ll _n0) {
    n0 = _n0;
    lim = sqrt(n0);
    ll fl = n0 / lim == lim;
    auto init = [&](ll n, int _id) {
        n %= M;
        sub(g1[_id] = n * (n + 1) / 2 % M, 1);
        sub(g2[_id] = n * (n + 1) % M * (n + n + 1) % M * i6 % M, 1);
    };
    auto doit = [&](ll n, int m, int i1, int i2) {
        sub(g1[i1], 1ll * (g1[i2] - s1[m - 1] + M) * p[m] % M);
        sub(g2[i1], 1ll * (g2[i2] - s2[m - 1] + M) * p2[m] % M);
    };
    For(i, 1, lim) init(i, i);
    For(i, 1, lim - fl) init(n0 / i, i + lim);
    For(i, 1, pc) {
        ll x = p[i], y = x * x, z = n0 / x;
        For(j, 1, min(lim - fl, n0 / y)) doit(n0 / j, i, j + lim, id(z / j));
        rFor(j, lim, min(lim + 1, y)) doit(j, i, j, j / x);
    }
}
int g(ll n) { return (g2[id(n)] - g1[id(n)] + M) % M; }
ll func(ll n) {
    n %= M;
    return n * (n - 1) % M;
}
int f(ll n, int m) {
    if(p[m] > n) return 0;
    ll ans = g(n) - s[m] + M;
    For(j, m + 1, pc) {
        int p = ::p[j];
        if(1ll * p * p > n) break;
        for(ll s = p; s * p <= n; s *= p) {
            ans += func(s) * f(n / s, j) % M + func(s * p);
        }
    }
    return ans % M;
}
ll solve(ll n) {
    init(), solveg(n);
    return (f(n, 0) + 1) % M;
}
