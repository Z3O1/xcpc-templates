// bostan_mori(): 多项式 · 线性递推第 n 项 (Bostan-Mori)

// 系数在 [0,M)，M=998244353；依赖 mul、ksm，O(k log k log(n+2))，k 为输入长度上界。
// k 阶递推：Q=1-Σ c_i x^i，P=(Q * Σ_{i<k} a_i x^i) mod x^k。

using poly = vector<int>;
// 返回 [x^n]P/Q；n>=0、Q 非空且 Q[0]!=0；P 可空或长于 Q。
int bostan_mori(poly P, poly Q, ll n) {
    if(P.empty()) return 0;
    if(Q.size() == 1) return n < (ll)P.size() ? 1ll * P[n] * ksm(Q[0], M - 2, M) % M : 0;
    while(n) {
        // 只求第 n 项,高于 n 次的分子/分母项可安全丢弃。
        if(n < (ll) P.size()) P.resize(n + 1);
        if(n < (ll) Q.size()) Q.resize(n + 1);
        poly Qm = Q;
        For(i, 1, (int) Qm.size() - 1) if(i & 1) Qm[i] = (M - Qm[i]) % M; // Q(-x)
        poly U = mul(P, Qm), V = mul(Q, Qm);
        P.clear();
        for(int i = n & 1; i < (int) U.size(); i += 2) P.push_back(U[i]); // 取偶数/奇数项
        Q.clear();
        for(int i = 0; i < (int) V.size(); i += 2) Q.push_back(V[i]);
        n >>= 1;
    }
    return 1ll * P[0] * ksm(Q[0], M - 2, M) % M;
}
