// rmq.cpp 自测:与逐段扫描取 min/max 对照(含重复值、单元素、块边界、dt=0/1 两种下标约定)
//
// 覆盖 rmq_t 的三个 bd 重载与 que:
//   bd(n, T*)        1-based 输入(dt=1)
//   bd(vector<T>&)   0-based 输入，查询下标仍从 dt 开始
//   bd(n, func)      func(i) 从偏移处取值
//   dt = 0 的 0-based 用法(rmq_t<int, less<int>, 0>)
#include "../_check_base.hpp"
#include "rmq.cpp"

// 覆盖三种构建入口与 dt=0/1 下标约定；失败返回非 0。
int main() {
    // 1) 1-based(dt=1):n = 1..80 与 2B±1 等块边界,全区间对穷举
    {
        ll cnt = 0;
        vector<int> sizes{1, 2, 3, 19, 20, 21, 39, 40, 41, 59, 60, 61, 100, 121, 200};
        for(int n : sizes) for(int mode = 0; mode < 4; ++mode) {
            vector<int> a(n + 1); // 1-based
            vector<int> raw(n);
            For(i, 1, n) {
                if(mode == 0) raw[i - 1] = (int) rnd(-1000, 1000);
                else if(mode == 1) raw[i - 1] = 7;                          // 全相等
                else if(mode == 2) raw[i - 1] = (int) rnd(0, 1);            // 大量重复
                else raw[i - 1] = (int) rnd(0, 5) * (i % 7 == 0 ? -1 : 1);  // 稀疏重复
                a[i] = raw[i - 1];
            }
            rmq_t<int> rmax;
            rmax.bd(n, a.data());               // 1-based 指针
            rmq_t<int, greater<int>> rmin;
            rmin.bd(n, a.data());
            For(l, 1, n) For(rr, l, n) {
                int mx = raw[l - 1], mn = raw[l - 1];
                For(i, l - 1, rr - 1) mx = max(mx, raw[i]), mn = min(mn, raw[i]);
                if(rmax.que(l, rr) != mx) {
                    printf("  [FAIL] max n=%d [%d,%d] want=%d got=%d\n", n, l, rr, mx, rmax.que(l, rr));
                    return 1;
                }
                if(rmin.que(l, rr) != mn) {
                    printf("  [FAIL] min n=%d [%d,%d] want=%d got=%d\n", n, l, rr, mn, rmin.que(l, rr));
                    return 1;
                }
                cnt += 2;
            }
        }
        printf("  [ok] 1-based dt=1:15 种规模 × 4 种数据 × 全部区间(%lld 次查询)\n", cnt);
    }

    // 2) bd(vector&) 重载(dt=1)
    {
        For(t, 1, 200) {
            int n = (int) rnd(1, 300);
            vector<int> v(n);
            For(i, 0, n - 1) v[i] = (int) rnd(-100, 100);
            rmq_t<int> rmax;
            rmax.bd(v);
            For(q, 1, 200) {
                int l = (int) rnd(1, n), rr = (int) rnd(l, n);
                int mx = v[l - 1];
                For(i, l - 1, rr - 1) mx = max(mx, v[i]);
                if(rmax.que(l, rr) != mx) return printf("  [FAIL] bd(vector) n=%d [%d,%d]\n", n, l, rr), 1;
            }
        }
        ok("200 组 bd(vector<int>&)(dt=1)一致");
    }

    // 3) dt = 0 的 0-based 用法
    {
        For(t, 1, 200) {
            int n = (int) rnd(1, 300);
            vector<int> v(n);
            For(i, 0, n - 1) v[i] = (int) rnd(-100, 100);
            rmq_t<int, less<int>, 0> rmax;
            rmax.bd(n, v.data());
            For(q, 1, 200) {
                int l = (int) rnd(0, n - 1), rr = (int) rnd(l, n - 1);
                int mx = v[l];
                For(i, l, rr) mx = max(mx, v[i]);
                if(rmax.que(l, rr) != mx) return printf("  [FAIL] dt=0 n=%d [%d,%d]\n", n, l, rr), 1;
            }
        }
        ok("200 组 dt=0 的 0-based 用法一致");
    }

    // 4) bd(n, func) 重载(dt=1:func 收到 1-based 下标)
    {
        For(t, 1, 100) {
            int n = (int) rnd(1, 200);
            rmq_t<int> rmax;
            rmax.bd(n, [&](int i) { return (int) ((i * 2654435761u) % 1000); });
            For(q, 1, 100) {
                int l = (int) rnd(1, n), rr = (int) rnd(l, n);
                int mx = INT_MIN;
                For(i, l, rr) mx = max(mx, (int) ((i * 2654435761u) % 1000));
                if(rmax.que(l, rr) != mx) return printf("  [FAIL] bd(n,func) n=%d [%d,%d]\n", n, l, rr), 1;
            }
        }
        ok("100 组 bd(n, func) 一致");
    }

    // 5) 单元素与全区间、单调数据的边界
    {
        vector<int> a1{0, 42}; // 1-based,只有 a[1]
        rmq_t<int> r;
        r.bd(1, a1.data());
        CHECK(r.que(1, 1) == 42, "n=1 单元素查询");
        vector<int> inc(1001), dec(1001);
        For(i, 1, 1000) inc[i] = i, dec[i] = 1000 - i + 1;
        rmq_t<int> ri, rd, rmin;
        ri.bd(1000, inc.data()), rd.bd(1000, dec.data());
        CHECK(ri.que(1, 1000) == 1000 && ri.que(500, 1000) == 1000 && ri.que(1, 500) == 500, "递增数组区间最大值");
        CHECK(rd.que(1, 1000) == 1000 && rd.que(2, 700) == 999, "递减数组区间最大值");
        CHECK(ri.que(1000, 1000) == 1000 && rd.que(1, 1) == 1000, "端点单点");
        rmq_t<int> big;
        big.bd(1000, inc.data());
        For(t, 1, 5000) {
            int l = (int) rnd(1, 1000), rr = (int) rnd(l, 1000);
            if(big.que(l, rr) != rr) return printf("  [FAIL] 递增数组 [%d,%d] 应为 %d\n", l, rr, rr), 1;
        }
        ok("单元素 / 递增 / 递减 / 端点边界一致");
    }

    // 任意正负 dt：三入口反复重建并逐区间对拍，回调必须恰好读取 n 次。
    auto offsets = [&]<int dt>() {
        rmq_t<int, less<int>, dt> r;
        for(int n : {1, 21, 80, 3, 41, 1}) {
            vector<int> v(n);
            for(int &x : v) x = rnd(-1000, 1000);
            auto verify = [&]() {
                For(l, 0, n - 1) For(rr, l, n - 1) {
                    int want = *max_element(v.begin() + l, v.begin() + rr + 1);
                    if(r.que(l + dt, rr + dt) != want) {
                        printf("[FAIL] rebuild dt=%d n=%d [%d,%d]\n", dt, n, l, rr);
                        for(int x : v) printf("%d ", x);
                        exit(1);
                    }
                }
            };
            r.bd(v), verify();
            int calls = 0;
            r.bd(n, [&](int i) {
                if(i != calls++ + dt) {
                    printf("[FAIL] callback dt=%d n=%d i=%d calls=%d\n", dt, n, i, calls);
                    exit(1);
                }
                return v[i - dt];
            });
            CHECK(calls == n, "回调恰好调用 n 次");
            verify();
            vector<int> padded(n + abs(dt));
            copy(all(v), padded.begin() + max(dt, 0));
            r.bd(n, padded.data() + max(-dt, 0)), verify();
        }
    };
    offsets.operator()<0>(), offsets.operator()<1>();
    offsets.operator()<7>(), offsets.operator()<-7>();
    PASSED("rmq");
}
