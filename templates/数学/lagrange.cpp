// lag()/lag_i()/lag_p(): 数学 · 拉格朗日插值

vector<mint> lag(vector<pair<mint, mint>> a) {
    int n = a.size();
    vector<mint> f(n + 1), ans(n + 1, 0);
    f[0] = 1;
    For(i, 1, n) {
        auto [x, y] = a[i - 1];
        rFor(j, i - 1, 0) f[j + 1] += f[j], f[j] *= -x;
    }
    For(i, 1, n) {
        auto [x, y] = a[i - 1];
        mint w = 1, inv = !x ? 0 : -x.inv(), z = 0;
        For(j, 1, n) if(i != j) w *= x - a[j - 1].first;
        w = w.inv() * y;
        if(!x) For(j, 1, n) ans[j - 1] += f[j] * w;
        else For(j, 0, n) z = (f[j] - z) * inv, ans[j] += z * w;
    }
    return ans;
}

// 需要 fac(n), ifac(n)
void lag_i(int n, mint *a) {
    vector<mint> f(n + 1, 0), g(n, 0);
    f[0] = 1;
    For(i, 1, n) rFor(j, i - 1, 0) f[j + 1] += f[j], f[j] *= mint::raw(mint::getM() - i);
    For(i, 1, n) {
        mint inv = -ifac(i) * fac(i - 1), w = a[i] * ifac(n - i) * ifac(i - 1) * ((i ^ n) & 1 ? -1 : 1), s = 0;
        For(j, 0, n - 1) s = (f[j] - s) * inv, g[j] += s * w;
    }
    For(i, 0, n - 1) a[i] = g[i];
    a[n] = 0;
}
mint lag_p(int n, mint *a, mint x) {
    if(x.val() < n) return a[x.val()];
    vector<mint> f(n + 1), g(n + 1);
    mint ans = 0;
    f[0] = g[n] = 1;
    For(i, 0, n - 1) f[i + 1] = f[i] * (x - i);
    rFor(i, n - 1, 0) g[i] = g[i + 1] * (x - i);
    For(i, 0, n - 1) {
        mint w = f[i] * g[i + 1] * a[i] * ifac(n - 1 - i) * ifac(i) * ((i ^ n) & 1 ? 1 : -1);
        ans += w;
    }
    return ans;
}
