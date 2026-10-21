// PN::solve(): 数论 · Powerful Number 筛

// 求积性函数前缀和 sum_{i=1}^n f(i) mod M，1 <= n <= 1e12；预处理与搜索均 O(sqrt(n))。
// 按题目更换 f_pe(p^e)、g_pe(p^e)（须满足 f(p)=g(p)）与 g 的前缀和 G(n)。
namespace PN {
const int N = 1e6 + 10, i2 = (M + 1) / 2;
int p[N / 10], pc, h[N / 10][42];
bool b[N];
ll n0, ans;
// 示例：f(p^e) = p^(2e-1)，素数处 f(p) = p = g(p)，伴随函数取 g(n) = n。
int f_pe(ll q, int e, ll qe) { return qe % M * qpow(q, e - 1) % M; }
int g_pe(ll q, int e, ll qe) { return qe % M; }
int G(ll n) { return n %= M, n * (n + 1) % M * i2 % M; }
void init() {
    if(pc) return;
    For(i, 2, N - 1) {
        if(!b[i]) p[++pc] = i;
        for(int j = 1, x; j <= pc && (x = i * p[j]) < N; ++j) {
            b[x] = 1;
            if(i % p[j] == 0) break;
        }
    }
}
void dfs(int k, ll d, ll hd) {
    ans = (ans + hd * G(n0 / d)) % M;
    For(i, k, pc) {
        ll q = p[i];
        if(q > n0 / q / d) break;
        ll x = d * q * q;
        for(int e = 2;; ++e, x *= q) {
            if(h[i][e]) dfs(i + 1, x, hd * h[i][e] % M);
            if(x > n0 / q) break;
        }
    }
}
ll solve(ll n) {
    init(), n0 = n, ans = 0;
    For(i, 1, pc) {
        ll q = p[i], qe = q;
        if(q * q > n) break;
        int gp[42];
        h[i][0] = 1, h[i][1] = 0, gp[1] = g_pe(q, 1, q);
        for(int e = 2; qe <= n / q; ++e) {
            qe *= q, gp[e] = g_pe(q, e, qe);
            ll v = f_pe(q, e, qe);
            For(j, 1, e) v = (v - 1ll * gp[j] * h[i][e - j]) % M;
            h[i][e] = (v + M) % M;
        }
    }
    return dfs(1, 1, 1), ans;
}
}  // namespace PN
