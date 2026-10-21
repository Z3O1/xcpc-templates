// 单纯形 自测:顶点枚举参考、LP强对偶定理2000组对拍、退化/Beale/Klee-Minty、无解与无界判定
#include "../_check_base.hpp"
#include "单纯形.cpp"

using namespace LP;

// 校验 x 是否满足 Ax <= b, x >= 0，且 c^T x == val。
static bool check_primal_feasible(const vector<vector<db>> &A, const vector<db> &b,
                                  const vector<db> &c, db val, const vector<db> &x) {
    int m = b.size(), n = c.size();
    if((int) x.size() != n) return false;
    For(j, 0, n - 1) if(x[j] < -1e-7L || !isfinite(x[j])) return false;
    For(i, 0, m - 1) {
        db s = 0;
        For(j, 0, n - 1) s += A[i][j] * x[j];
        if(s > b[i] + 1e-7L * (1 + abs(b[i]))) return false;
    }
    db obj = 0;
    For(j, 0, n - 1) obj += c[j] * x[j];
    return abs(obj - val) <= 1e-7L * (1 + abs(val));
}

// 求解对偶问题: min b^T y s.t. A^T y >= c, y >= 0
// 等价于: - max (-b)^T y s.t. (-A^T) y <= -c, y >= 0
static tuple<int, db, vector<db>> solve_dual(const vector<vector<db>> &A, const vector<db> &b, const vector<db> &c) {
    int m = b.size(), n = c.size();
    vector<vector<db>> AT(n, vector<db>(m));
    vector<db> nb(m), nc(n);
    For(i, 0, m - 1) nb[i] = -b[i];
    For(j, 0, n - 1) {
        nc[j] = -c[j];
        For(i, 0, m - 1) AT[j][i] = -A[i][j];
    }
    auto [st, v, y] = simplex(AT, nc, nb);
    return {st, -v, y};
}

// 2D/3D 小规模精确 Cramer 法枚举所有基本解，作为独立参考。
static pair<int, db> brute_vertices_2d(const vector<vector<db>> &A, const vector<db> &b, const vector<db> &c) {
    int m = b.size();
    // 收集所有边界线: A[i][0] x + A[i][1] y = b[i], 以及 x = 0, y = 0
    vector<array<db, 3>> lines;
    For(i, 0, m - 1) lines.push_back({A[i][0], A[i][1], b[i]});
    lines.push_back({-1, 0, 0});
    lines.push_back({0, -1, 0});
    int L = lines.size();
    bool found = false;
    db best = -1e100L;
    For(i, 0, L - 1) For(j, i + 1, L - 1) {
        db det = lines[i][0] * lines[j][1] - lines[i][1] * lines[j][0];
        if(abs(det) < 1e-12L) continue;
        db x = (lines[i][2] * lines[j][1] - lines[i][1] * lines[j][2]) / det;
        db y = (lines[i][0] * lines[j][2] - lines[i][2] * lines[j][0]) / det;
        if(x < -1e-9L || y < -1e-9L) continue;
        bool ok = true;
        For(k, 0, m - 1) if(A[k][0] * x + A[k][1] * y > b[k] + 1e-9L) { ok = false; break; }
        if(ok) {
            found = true;
            cmax(best, c[0] * x + c[1] * y);
        }
    }
    return {found ? 1 : 0, best};
}

int main() {
    // 1. 解析用例：经典有界、b<0 需 Phase 1、无解、无界、空维度
    {
        // max 3x + 4y s.t. x + y <= 4, x + 3y <= 6, x,y >= 0 -> 在 (3, 1) 取 13
        auto [st, val, x] = simplex({{1, 1}, {1, 3}}, {4, 6}, {3, 4});
        CHECK(st == 1 && abs(val - 13) < 1e-9L && abs(x[0] - 3) < 1e-9L && abs(x[1] - 1) < 1e-9L, "2D 标准有界最优");

        // 含 b_i < 0 的可行问题: x + y >= 2 (-x - y <= -2), x <= 3, y <= 3, max x + 2y -> (3, 3) 取 9
        auto [st2, val2, x2] = simplex({{-1, -1}, {1, 0}, {0, 1}}, {-2, 3, 3}, {1, 2});
        CHECK(st2 == 1 && abs(val2 - 9) < 1e-9L && check_primal_feasible({{-1, -1}, {1, 0}, {0, 1}}, {-2, 3, 3}, {1, 2}, val2, x2), "b<0 Phase 1 可行");

        // 无可行解: x + y <= 1 且 x + y >= 2 (-x - y <= -2)
        auto [st3, val3, x3] = simplex({{1, 1}, {-1, -1}}, {1, -2}, {1, 1});
        CHECK(st3 == 0 && x3.empty(), "无可行解判定");

        // 无界: max x + y s.t. x - y <= 1, -x + y <= 1
        auto [st4, val4, x4] = simplex({{1, -1}, {-1, 1}}, {1, 1}, {1, 1});
        CHECK(st4 == -1 && val4 == inf && x4.empty(), "无界判定");

        // 空矩阵或 0 变量边界
        auto [st5, val5, x5] = simplex({}, {}, {1, 2});
        CHECK(st5 == -1 && val5 == inf && x5.empty(), "无约束正目标无界");
        auto [st6, val6, x6] = simplex({}, {}, {-1, -2});
        CHECK(st6 == 1 && abs(val6) < 1e-9L && x6 == vector<db>({0, 0}), "无约束负目标在原点最优");
        auto [st7, val7, x7] = simplex({{}, {}}, {1, 0}, {});
        CHECK(st7 == 1 && val7 == 0 && x7.empty(), "n=0 可行");
        auto [st8, val8, x8] = simplex({{}, {}}, {1, -1}, {});
        CHECK(st8 == 0, "n=0 不可行");
    }
    ok("解析用例：有界/Phase 1/无解/无界/0维边界");

    // 2. 经典退化防循环测试 (Beale 退化反例 + Klee-Minty 立方体)
    {
        // Beale 1955 经典退化例子:
        // max 0.75 x1 - 150 x2 + 0.02 x3 - 6 x4
        // s.t. 0.25 x1 - 60 x2 - 0.04 x3 + 9 x4 <= 0
        //      0.50 x1 - 90 x2 - 0.02 x3 + 3 x4 <= 0
        //      x3 <= 1
        vector<vector<db>> A = {
            {0.25L, -60.0L, -0.04L, 9.0L},
            {0.50L, -90.0L, -0.02L, 3.0L},
            {0.0L, 0.0L, 1.0L, 0.0L},
        };
        vector<db> b = {0, 0, 1}, c = {0.75L, -150.0L, 0.02L, -6.0L};
        auto [st, val, x] = simplex(A, b, c);
        CHECK(st == 1 && abs(val - 0.05L) < 1e-9L && check_primal_feasible(A, b, c, val, x), "Beale 退化例子不死循环且精确到 1/20");

        // Klee-Minty 立方体 (d=6)
        int d = 6;
        vector<vector<db>> KA(d, vector<db>(d, 0));
        vector<db> Kb(d), Kc(d);
        For(i, 0, d - 1) {
            KA[i][i] = 1;
            For(j, 0, i - 1) KA[i][j] = ldexp(db(1), i - j + 1);
            Kb[i] = powl(5, i + 1);
            Kc[i] = ldexp(db(1), d - 1 - i);
        }
        auto [kst, kval, kx] = simplex(KA, Kb, Kc);
        CHECK(kst == 1 && abs(kval - powl(5, d)) < 1e-7L && check_primal_feasible(KA, Kb, Kc, kval, kx), "Klee-Minty 立方体");
    }
    ok("Beale 退化反例与 Klee-Minty 立方体");

    // 3. 800 组 2D 随机半平面（含包围盒保证有界，允许无解）与 Cramer 顶点枚举对拍
    For(t, 1, 800) {
        int m = rnd(2, 10);
        vector<vector<db>> A(m + 2, vector<db>(2));
        vector<db> b(m + 2), c = {db(rnd(-20, 20)), db(rnd(-20, 20))};
        For(i, 0, m - 1) {
            A[i][0] = rnd(-10, 10), A[i][1] = rnd(-10, 10);
            b[i] = rnd(-5, 30);
        }
        A[m] = {1, 0}, b[m] = 50;
        A[m + 1] = {0, 1}, b[m + 1] = 50;
        auto [st, val, x] = simplex(A, b, c);
        auto [wst, wval] = brute_vertices_2d(A, b, c);
        if(st != wst || (st == 1 && (!check_primal_feasible(A, b, c, val, x) || abs(val - wval) > 1e-7L * (1 + abs(wval))))) {
            printf("  [FAIL] 2D 顶点枚举对拍失败 case=%d st=%d wst=%d val=%Lg wval=%Lg\n", t, st, wst, val, wval);
            return 1;
        }
    }
    ok("800 组 2D 随机线性规划（含 b<0 与无解）与独立顶点枚举完全一致");

    // 4. 1200 组高维随机 LP 强对偶性检验 (原始可行有界 <=> 对偶可行有界 且 c^T x* == b^T y*)
    For(t, 1, 1200) {
        int n = rnd(1, 12), m = rnd(1, 15);
        vector<vector<db>> A(m, vector<db>(n));
        vector<db> b(m), c(n);
        For(i, 0, m - 1) {
            For(j, 0, n - 1) A[i][j] = rnd(-10, 10);
            b[i] = rnd(t % 3 == 0 ? -10 : 0, 30);
        }
        For(j, 0, n - 1) c[j] = rnd(-10, 10);
        auto [pst, pval, x] = simplex(A, b, c);
        auto [dst, dval, y] = solve_dual(A, b, c);
        if(pst == 1) {
            if(dst != 1 || !check_primal_feasible(A, b, c, pval, x) || abs(pval - dval) > 1e-6L * (1 + abs(pval))) {
                printf("  [FAIL] 强对偶不符 case=%d pst=%d dst=%d pval=%.15Lg dval=%.15Lg\n", t, pst, dst, pval, dval);
                return 1;
            }
        } else if(pst == -1) {
            if(dst != 0) {
                printf("  [FAIL] 原始无界但对偶非无解 case=%d dst=%d\n", t, dst);
                return 1;
            }
        } else {
            if(dst == 1) {
                printf("  [FAIL] 原始无解但对偶有有界解 case=%d\n", t);
                return 1;
            }
        }
    }
    ok("1200 组高维随机 LP (含 b<0、无解、无界) 原始-对偶强对偶定理校验通过");

    PASSED("单纯形");
}
