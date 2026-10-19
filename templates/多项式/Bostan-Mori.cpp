// bostan_mori(): 多项式 · 线性递推第 n 项 (Bostan-Mori)

// 依赖 6.3 NTT（mul）。
// 求 P/Q 的第 n 项系数 [x^n]P/Q，O(k log k log n)，k 为输入长度上界；系数在 [0,M)。

using poly = vector<int>;
// 返回 [x^n]P/Q；n>=0、Q 非空且 Q[0]!=0；P 可空或长于 Q。
int bostan_mori(poly p, poly q, ll n) {
    if(p.empty()) return 0;
    if(q.size() == 1) return n < (ll)p.size() ? 1ll * p[n] * qpow(q[0], M - 2) % M : 0;
    while(n) {
        if(n < (ll) p.size()) p.resize(n + 1);
        if(n < (ll) q.size()) q.resize(n + 1);
        poly Qm = q;
        For(i, 1, (int) Qm.size() - 1) if(i & 1) Qm[i] = (M - Qm[i]) % M;
        poly U = mul(p, Qm), V = mul(q, Qm);
        p.clear();
        for(int i = n & 1; i < (int) U.size(); i += 2) p.push_back(U[i]);
        q.clear();
        for(int i = 0; i < (int) V.size(); i += 2) q.push_back(V[i]);
        n >>= 1;
    }
    return 1ll * p[0] * qpow(q[0], M - 2) % M;
}
