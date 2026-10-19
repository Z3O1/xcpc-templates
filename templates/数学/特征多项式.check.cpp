// 特征多项式 自测:Leibniz 展开、Jordan/重根/零矩阵、随机对拍、点值行列式与 Cayley-Hamilton
#include "../_check_base.hpp"
#include "特征多项式.cpp"

static void fail(const mat &a, const char *what, const poly &got = {}, const poly &want = {}) {
    cerr << what << " n=" << a.size() << '\n';
    for(const auto &row : a) {
        for(int x : row) cerr << x << ' ';
        cerr << '\n';
    }
    cerr << "got: ";
    for(int x : got) cerr << x << ' ';
    cerr << "\nwant: ";
    for(int x : want) cerr << x << ' ';
    cerr << '\n';
    exit(1);
}

// 独立按排列定义 det(xI-a)，每项直接乘 n 个一次因子。
static poly ref_charpoly(const mat &a) {
    int n = a.size();
    vector<int> p(n);
    iota(all(p), 0);
    poly ans(n + 1);
    do {
        int parity = 0;
        For(i, 0, n - 1) For(j, 0, i - 1) parity ^= p[j] > p[i];
        poly term{1};
        For(i, 0, n - 1) {
            poly next(term.size() + 1);
            For(k, 0, (int)term.size() - 1) {
                next[k] = (next[k] + (ll)term[k] * (M - a[i][p[i]])) % M;
                if(p[i] == i) next[k + 1] = (next[k + 1] + term[k]) % M;
            }
            term = next;
        }
        For(k, 0, n) ans[k] = (ans[k] + (parity ? M - term[k] : term[k])) % M;
    } while(next_permutation(all(p)));
    return ans;
}

static int ref_det(mat a) {
    int n = a.size();
    ll ans = 1;
    For(i, 0, n - 1) {
        int pivot = i;
        while(pivot < n && !a[pivot][i]) ++pivot;
        if(pivot == n) return 0;
        if(pivot != i) swap(a[pivot], a[i]), ans = (M - ans) % M;
        ans = ans * a[i][i] % M;
        ll inv = qpow(a[i][i], M - 2);
        For(j, i + 1, n - 1) {
            ll c = a[j][i] * inv % M;
            For(k, i, n - 1) a[j][k] = (a[j][k] - c * a[i][k] % M + M) % M;
        }
    }
    return ans;
}

static mat multiply(const mat &a, const mat &b) {
    int n = a.size();
    mat c(n, poly(n));
    For(i, 0, n - 1) For(j, 0, n - 1) For(k, 0, n - 1)
        c[i][j] = (c[i][j] + (ll)a[i][k] * b[k][j]) % M;
    return c;
}

static void verify(const mat &a, bool leibniz, bool cayley) {
    int n = a.size();
    poly got = charpoly(a);
    if((int)got.size() != n + 1 || got[n] != 1) fail(a, "size/monic", got);
    for(int x : got)
        if(x < 0 || x >= M) fail(a, "coefficient range", got);
    if(leibniz) {
        poly want = ref_charpoly(a);
        if(got != want) fail(a, "Leibniz", got, want);
    }
    for(int x : {0, 1, M - 1, (int)rnd(0, M - 1)}) {
        mat b = a;
        For(i, 0, n - 1) For(j, 0, n - 1) b[i][j] = ((i == j ? x : 0) - a[i][j] + M) % M;
        int value = 0;
        rFor(i, n, 0) value = ((ll)value * x + got[i]) % M;
        int want = ref_det(b);
        if(value != want) {
            cerr << "x=" << x << '\n';
            fail(a, "point determinant", {value}, {want});
        }
    }
    if(cayley) {
        mat value(n, poly(n));
        rFor(k, n, 0) {
            value = multiply(value, a);
            For(i, 0, n - 1) value[i][i] = (value[i][i] + got[k]) % M;
        }
        for(const auto &row : value)
            for(int x : row)
                if(x) fail(a, "Cayley-Hamilton", got);
    }
    if(charpoly(a) != got) fail(a, "repeat call", got);
}

int main() {
    rng.seed(0x43484152504f4c59ULL);
    verify({}, true, true);
    for(int n : {1, 2, 3, 6, 12, 60}) {
        verify(mat(n, poly(n)), n <= 6, n <= 12);
        for(int root : {0, 1, M - 1, 314159}) {
            mat a(n, poly(n));
            For(i, 0, n - 1) {
                a[i][i] = root;
                if(i + 1 < n) a[i][i + 1] = 1;
            }
            // Jordan 块的独立闭式为 (x-root)^n。
            poly want{1};
            For(i, 1, n) {
                poly next(want.size() + 1);
                For(k, 0, (int)want.size() - 1) {
                    next[k] = (next[k] + (ll)want[k] * (M - root)) % M;
                    next[k + 1] = (next[k + 1] + want[k]) % M;
                }
                want = next;
            }
            if(charpoly(a) != want) fail(a, "Jordan closed form", charpoly(a), want);
            verify(a, n <= 6, n <= 12);
            For(i, 0, n - 1) For(j, 0, i - 1) swap(a[i][j], a[j][i]);
            verify(a, n <= 6, n <= 12);
        }
    }
    For(test, 0, 2499) {
        int n = rnd(0, 6);
        mat a(n, poly(n));
        For(i, 0, n - 1) For(j, 0, n - 1) {
            a[i][j] = test % 3 ? rnd(0, M - 1) : rnd(0, 2);
            if(test % 7 == 0 && i > j) a[i][j] = 0;
            if(test % 7 == 1 && rnd(0, 2)) a[i][j] = 0;
        }
        if(n > 1 && test % 7 == 2) a[n - 1] = a[0];
        verify(a, true, true);
    }
    For(test, 1, 60) {
        int n = test <= 40 ? 12 : 60;
        mat a(n, poly(n));
        For(i, 0, n - 1) For(j, 0, n - 1) a[i][j] = rnd(0, M - 1);
        verify(a, false, n == 12);
    }
    PASSED("特征多项式");
}
