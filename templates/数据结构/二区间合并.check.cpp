// 二区间合并 自测:非交换字符串/矩阵半群独立扫描、块边界/重建/无默认构造、合并次数界
#include "../_check_base.hpp"
#include "二区间合并.cpp"

struct test_value {
    ll x;
    explicit test_value(ll v) : x(v) {}
};

static array<int, 4> matrix_product(const array<int, 4> &a, const array<int, 4> &b) {
    array<int, 4> c{};
    For(i, 0, 1) For(j, 0, 1) For(k, 0, 1)
        c[2 * i + j] = (c[2 * i + j] + (ll)a[2 * i + k] * b[2 * k + j]) % M;
    return c;
}

int main() {
    rng.seed(20261019);
    auto concat = [](const string &a, const string &b) { return a + b; };
    disjoint_sparse_t<string, decltype(concat)> words({}, concat);
    vector<int> sizes{0, 1, 2, 3, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129};
    For(rep, 1, 300) {
        int n = rep <= (int)sizes.size() ? sizes[rep - 1] : rnd(1, 130);
        vector<string> a(n);
        for(auto &s : a) s = char('a' + rnd(0, 25));
        words.build(a);
        For(l, 0, n - 1) {
            string want;
            For(r, l, n - 1) {
                want += a[r];
                string got = words.query(l, r);
                if(got != want) {
                    printf("FAIL concat rep=%d n=%d [%d,%d] input=", rep, n, l, r);
                    for(auto s : a) printf("%s", s.c_str());
                    printf(" want=%s got=%s\n", want.c_str(), got.c_str());
                    return 1;
                }
            }
        }
    }
    ok("字符串全部非空区间：非交换、非幂等、二次幂边界及空表重建");

    For(rep, 1, 400) {
        int n = rnd(1, 300);
        vector<array<int, 4>> a(n);
        for(auto &p : a)
            for(int &x : p) x = rnd(0, M - 1);
        disjoint_sparse_t<array<int, 4>, decltype(&matrix_product)> table(a, matrix_product);
        For(q, 1, 200) {
            int l = rnd(0, n - 1);
            int r = rnd(l, n - 1);
            array<int, 4> want{1, 0, 0, 1};
            For(i, l, r) {
                auto b = want;
                want = {(int)(((ll)b[0] * a[i][0] + (ll)b[1] * a[i][2]) % M),
                        (int)(((ll)b[0] * a[i][1] + (ll)b[1] * a[i][3]) % M),
                        (int)(((ll)b[2] * a[i][0] + (ll)b[3] * a[i][2]) % M),
                        (int)(((ll)b[2] * a[i][1] + (ll)b[3] * a[i][3]) % M)};
            }
            auto got = table.query(l, r);
            if(got != want) {
                printf("FAIL matrix rep=%d n=%d [%d,%d]\n", rep, n, l, r);
                for(auto p : a) printf("%d %d %d %d\n", p[0], p[1], p[2], p[3]);
                for(int x : want) printf("want %d\n", x);
                for(int x : got) printf("got %d\n", x);
                return 1;
            }
        }
    }
    ok("400 组非交换矩阵半群，80000 次独立顺序扫描");

    ll calls = 0;
    auto sum = [&](const test_value &a, const test_value &b) {
        ++calls;
        return test_value(a.x + b.x);
    };
    int n = 200003;
    vector<test_value> a;
    For(i, 0, n - 1) a.emplace_back(i + 1);
    disjoint_sparse_t<test_value, decltype(sum)> table(a, sum);
    CHECK(calls <= 2LL * n * (32 - __builtin_clz(n)), "预处理合并次数 O(n log n)，无需默认构造");
    For(q, 1, 200000) {
        int l = rnd(0, n - 1);
        int r = rnd(l, n - 1);
        calls = 0;
        ll got = table.query(l, r).x;
        ll want = (ll)(l + r + 2) * (r - l + 1) / 2;
        if(got != want || calls != (l != r)) {
            printf("FAIL sum [%d,%d] want=%lld got=%lld calls=%lld\n", l, r, want, got, calls);
            return 1;
        }
    }
    ok("20 万规模与 20 万次查询：单点零合并，其余恰好一次");
    PASSED("二区间合并");
}
