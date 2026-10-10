// Berlekamp-Massey 自测:已知递推还原、独立扩欧+高斯消元求最短阶数、穷举短序列、退化及长序列
#include "../_check_base.hpp"
#include "Berlekamp-Massey.cpp"

// 独立参考逆元:扩欧,不复用模板的快速幂。
static int ref_inv(int v) {
    ll a = v, b = MOD, x = 1, y = 0;
    while(b) {
        ll q = a / b;
        swap(a -= q * b, b), swap(x -= q * y, y);
    }
    return (x % MOD + MOD) % MOD;
}
struct LS {
    bool cons, uniq;
    poly sol;
};
static LS brute_order(const poly &a, int k) {
    vector<poly> A;
    poly b;
    For(i, k, (int)a.size() - 1) {
        poly row(k);
        For(j, 1, k) row[j - 1] = a[i - j];
        A.push_back(row), b.push_back(a[i]);
    }
    int r = A.size(), row = 0;
    vector<int> where(k, -1);
    for(int col = 0; col < k && row < r; ++col) {
        int piv = row;
        while(piv < r && !A[piv][col]) ++piv;
        if(piv == r) continue;
        swap(A[row], A[piv]), swap(b[row], b[piv]);
        int iv = ref_inv(A[row][col]);
        For(j, 0, k - 1) A[row][j] = (ll)A[row][j] * iv % MOD;
        b[row] = (ll)b[row] * iv % MOD;
        For(i, 0, r - 1) if(i != row && A[i][col]) {
            int f = A[i][col];
            For(j, 0, k - 1) A[i][j] = (A[i][j] - (ll)f * A[row][j] % MOD + MOD) % MOD;
            b[i] = (b[i] - (ll)f * b[row] % MOD + MOD) % MOD;
        }
        where[col] = row++;
    }
    For(i, row, r - 1) if(b[i]) return {false, false, {}};
    poly sol(k);
    For(i, 0, k - 1) if(where[i] >= 0) sol[i] = b[where[i]];
    return {true, row == k, sol};
}
static string pstr(const poly &a) {
    string s = "[" + to_string(a.size()) + "]";
    For(i, 0, min(20, (int)a.size()) - 1) s += " " + to_string(a[i]);
    return s;
}
static void fail(const char *what, const poly &a, const poly &got, const poly &want = {}) {
    cerr << "[FAIL] " << what << "\na=" << pstr(a) << "\ngot=" << pstr(got)
         << "\nwant=" << pstr(want) << endl;
    exit(1);
}
static void fits(const poly &a, const poly &c) {
    For(i, 0, (int)c.size() - 1) if(c[i] < 0 || c[i] >= MOD) fail("系数未归一", a, c);
    if(c.size() > a.size()) fail("阶数超过项数", a, c);
    For(i, (int)c.size(), (int)a.size() - 1) {
        ll s = 0;
        For(j, 1, (int)c.size()) s = (s + (ll)c[j - 1] * a[i - j]) % MOD;
        if(s != a[i]) fail("递推不能拟合", a, c);
    }
}
static void verify(const poly &a) {
    poly got = BM(a);
    fits(a, got);
    For(k, 0, (int)a.size()) {
        LS ref = brute_order(a, k);
        if(!ref.cons) continue;
        if((int)got.size() != k) fail("不是最短阶数", a, got, ref.sol);
        if(ref.uniq && got != ref.sol) fail("唯一解不同", a, got, ref.sol);
        return;
    }
    fail("参考无解", a, got);
}
static poly random_poly(int n) {
    poly a(n);
    for(int &x : a) x = rnd(0, MOD - 1);
    return a;
}
static poly gen(poly a, const poly &c, int n) {
    int k = c.size();
    a.resize(max(n, k));
    For(i, k, n - 1) {
        ll s = 0;
        For(j, 1, k) s = (s + (ll)c[j - 1] * a[i - j]) % MOD;
        a[i] = s;
    }
    a.resize(n);
    return a;
}
int main() {
    rng.seed(20040924);
    For(t, 1, 3000) {
        int k = rnd(1, 12);
        poly c = random_poly(k), a0 = random_poly(k);
        c.back() = rnd(1, MOD - 1), a0[0] = rnd(1, MOD - 1);
        poly a = gen(a0, c, 2 * k + rnd(0, 4)), got = BM(a);
        fits(a, got);
        if(got != c) fail("随机递推未还原", a, got, c);
    }
    For(t, 1, 200) {
        int k = rnd(13, 40);
        poly c(k), a0 = random_poly(k);
        for(int &x : c) x = rnd(0, 2);
        c.back() = 1;
        poly a = gen(a0, c, 2 * k + rnd(0, 10)), got = BM(a);
        if(got != c) fail("稀疏递推未还原", a, got, c);
    }
    ok("3000 组随机递推 + 200 组稀疏递推还原");
    For(t, 1, 1500) {
        poly a = random_poly(rnd(1, 16));
        if(t % 5 == 0) for(int &x : a) x = rnd(0, 1);
        verify(a);
    }
    For(t, 1, 500) {
        int k = rnd(2, 10), n = rnd(1, 2 * k - 1);
        poly a = gen(random_poly(k), random_poly(k), n);
        fits(a, BM(a));
        verify(a);
    }
    // 全枚举保证覆盖 MOD-1 的加减乘法,不只靠小数和随机。
    For(n, 0, 7) {
        int lim = 1;
        For(i, 1, n) lim *= 3;
        For(mask, 0, lim - 1) {
            poly a(n);
            int x = mask;
            for(int &v : a) v = (x % 3 == 2 ? MOD - 1 : x % 3), x /= 3;
            verify(a);
        }
    }
    ok("1500 组随机短序列 + 500 组不足前缀 + 长度≤7 的三值穷举与高斯消元一致");
    For(n, 0, 30) {
        if(!BM(poly(n)).empty()) fail("全零序列", poly(n), BM(poly(n)));
        if(n >= 2) for(int x : {1, MOD - 1}) {
            poly a(n, x), got = BM(a);
            if(got != poly{1}) fail("常数列", a, got, {1});
        }
    }
    For(t, 1, 20) {
        int r = rnd(1, MOD - 1);
        poly a = gen({1}, {r}, 20), got = BM(a);
        if(got != poly{r}) fail("等比列", a, got, {r});
    }
    for(auto [a0, c] : vector<pair<poly, poly>>{{{0, 1}, {1, 1}}, {{1, 2, 3}, {1, 0, 1}}}) {
        poly a = gen(a0, c, 30), got = BM(a);
        if(got != c) fail("经典递推", a, got, c);
    }
    For(t, 1, 200) {
        poly a(rnd(1, 20));
        For(i, 0, min((int)a.size(), (int)rnd(1, 5)) - 1) a[i] = rnd(1, MOD - 1);
        verify(a);
    }
    For(n, 1, 20) for(int x : {1, MOD - 1}) {
        poly a(n);
        a.back() = x;
        if(BM(a).size() != a.size()) fail("末项非零的阶数/越界回归", a, BM(a));
        verify(a);
    }
    for(int x : {1, MOD - 1}) if(BM({x}) != poly{0}) fail("单项非零的返回约定", {x}, BM({x}), {0});
    ok("空/零/常数/等比/经典递推/幂零/末项非零/单项边界");
    poly c = random_poly(60), a0 = random_poly(60);
    c.back() = rnd(1, MOD - 1), a0[0] = 1;
    for(int n : {200, 2000}) {
        poly a = gen(a0, c, n), got = BM(a);
        if(got != c) fail("长序列", a, got, c);
    }
    PASSED("Berlekamp-Massey");
}
