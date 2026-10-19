// sieve_t: 数论 · 线性筛 (素数 + 欧拉函数 + 莫比乌斯函数 + 约数个数 + 约数和)
// hide

// 依赖 全局 N。
// 线性筛 [1,N]，同时算出欧拉函数、莫比乌斯函数、约数个数与约数和；O(N)。
struct sieve_t {
    int p[N / 10 + 5], tot, mn[N + 1], phi[N + 1], mu[N + 1], d[N + 1], cnt[N + 1], lpf[N + 1];
    ll s[N + 1];
    // 在初始全零对象上筛 [1,N]，同一对象仅调用一次；要求 N>=1、2N<=INT_MAX。
    // 素数表 1-based，须满足 π(N) <= N/10+4。
    void work() {
        tot = 0;
        phi[1] = mu[1] = d[1] = s[1] = 1;
        For(i, 2, N) {
            if(!mn[i]) {
                p[++tot] = i;
                mn[i] = i;
                phi[i] = i - 1;
                mu[i] = -1;
                d[i] = 2;
                cnt[i] = 1;
                s[i] = i + 1;
                lpf[i] = i + 1;
            }
            for(int j = 1; j <= tot; ++j) {
                int x = i * p[j];
                if(x > N) break;
                mn[x] = p[j];
                cnt[x] = 1;
                lpf[x] = p[j] + 1;
                if(i % p[j] == 0) {
                    phi[x] = phi[i] * p[j];
                    mu[x] = 0;
                    cnt[x] = cnt[i] + 1;
                    lpf[x] = lpf[i] * p[j] + 1;
                    d[x] = d[i] / (cnt[i] + 1) * (cnt[x] + 1);
                    s[x] = s[i] / lpf[i] * lpf[x];
                    break;
                }
                phi[x] = phi[i] * (p[j] - 1);
                mu[x] = -mu[i];
                d[x] = d[i] * 2;
                s[x] = s[i] * (p[j] + 1);
            }
        }
    }
    // 要求 work() 已完成且 0 <= x <= N；返回 x 是否为素数。
    bool isp(int x) const { return x >= 2 && mn[x] == x; }
} S;
// 莫比乌斯/欧拉函数的前缀和(杜教筛要用):筛完后
//   For(i, 1, N) smu[i] = smu[i - 1] + S.mu[i], sphi[i] = sphi[i - 1] + S.phi[i];
