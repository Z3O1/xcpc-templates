// 行列式多项式 自测:独立 Leibniz 展开、双奇异/降次/零矩阵、随机对拍、点值与系数反转
#include "../_check_base.hpp"
#include "特征多项式.cpp"
#include "行列式多项式.cpp"

static void fail(const mat &a, const mat &b, const char *what, const poly &got = {}, const poly &want = {}) {
    cerr << what << " n=" << a.size() << '\n';
    cerr << "a:\n";
    for(const auto &row : a) {
        for(int x : row) cerr << x << ' ';
        cerr << '\n';
    }
    cerr << "b:\n";
    for(const auto &row : b) {
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

// 按排列定义逐项相乘一次多项式，不调用 charpoly 或待测消元。
static poly ref_det_poly(const mat &a, const mat &b) {
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
                next[k] = (next[k] + (ll)term[k] * b[i][p[i]]) % M;
                next[k + 1] = (next[k + 1] + (ll)term[k] * a[i][p[i]]) % M;
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
        ll inv = qpow(a[i][i], M - 2, M);
        For(j, i + 1, n - 1) {
            ll c = a[j][i] * inv % M;
            For(k, i, n - 1) a[j][k] = (a[j][k] - c * a[i][k] % M + M) % M;
        }
    }
    return ans;
}

static void verify(const mat &a, const mat &b, bool leibniz) {
    int n = a.size();
    poly got = det_poly(a, b);
    if((int)got.size() != n + 1) fail(a, b, "size", got);
    for(int x : got)
        if(x < 0 || x >= M) fail(a, b, "coefficient range", got);
    if(leibniz) {
        poly want = ref_det_poly(a, b);
        if(got != want) fail(a, b, "Leibniz", got, want);
    }
    if(got[0] != ref_det(b) || got[n] != ref_det(a)) fail(a, b, "end coefficients", got);
    poly reversed = det_poly(b, a);
    reverse(all(reversed));
    if(got != reversed) fail(a, b, "reverse coefficients", got, reversed);
    for(int x : {0, 1, M - 1, (int)rnd(0, M - 1)}) {
        mat c = b;
        For(i, 0, n - 1) For(j, 0, n - 1) c[i][j] = ((ll)a[i][j] * x + b[i][j]) % M;
        int value = 0;
        rFor(i, n, 0) value = ((ll)value * x + got[i]) % M;
        int want = ref_det(c);
        if(value != want) {
            cerr << "x=" << x << '\n';
            fail(a, b, "point determinant", {value}, {want});
        }
    }
    if(det_poly(a, b) != got) fail(a, b, "repeat call", got);
}

int main() {
    rng.seed(0x444554504f4c59ULL);
    verify({}, {}, true);
    for(int n : {1, 2, 3, 6, 12, 60}) {
        mat zero(n, poly(n));
        mat id = zero;
        mat shift = zero;
        For(i, 0, n - 1) {
            id[i][i] = 1;
            if(i + 1 < n) shift[i][i + 1] = 1;
        }
        verify(zero, zero, n <= 6);
        verify(zero, id, n <= 6);
        verify(id, zero, n <= 6);
        verify(shift, id, n <= 6);  // det(I+xN)=1，连续移次达到 n 次。
        verify(id, shift, n <= 6);
        verify(shift, shift, n <= 6);
        verify(id, id, n <= 6);
        // n>=2 时 a,b 均奇异，但 det(ax+b)=x*(x+1)^(n-2) 并非零多项式。
        mat a = id;
        mat b = id;
        a[0][0] = 0;
        if(n > 1) b[1][1] = 0;
        verify(a, b, n <= 6);
        // 列交换枢轴与重根。
        reverse(all(a));
        verify(a, b, n <= 6);
    }
    For(test, 0, 2499) {
        int n = rnd(0, 6);
        mat a(n, poly(n));
        mat b = a;
        For(i, 0, n - 1) For(j, 0, n - 1) {
            a[i][j] = test % 3 ? rnd(0, M - 1) : rnd(0, 2);
            b[i][j] = test % 3 ? rnd(0, M - 1) : rnd(0, 2);
            if(test % 7 == 0 && i > j) a[i][j] = b[i][j] = 0;
            if(test % 7 == 1 && rnd(0, 2)) a[i][j] = 0;
        }
        if(n && test % 7 == 2) fill(all(a[n - 1]), 0);
        if(n && test % 7 == 3) fill(all(b[n - 1]), 0);
        if(n > 1 && test % 7 == 4) a[n - 1] = a[0], b[n - 1] = b[0];
        verify(a, b, true);
    }
    For(test, 1, 60) {
        int n = test <= 40 ? 12 : 60;
        mat a(n, poly(n));
        mat b = a;
        For(i, 0, n - 1) For(j, 0, n - 1) a[i][j] = rnd(0, M - 1), b[i][j] = rnd(0, M - 1);
        if(test % 2 == 0) fill(all(a[n - 1]), 0);
        if(test % 3 == 0) fill(all(b[0]), 0);
        verify(a, b, false);
    }
    PASSED("行列式多项式");
}
