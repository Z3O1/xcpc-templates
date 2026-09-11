// bostan_mori(): 多项式 · 线性递推第 n 项 (Bostan-Mori)
// 求 [x^n] P(x) / Q(x)(P/Q 展开成幂级数后的第 n 项系数)。要求 Q[0] != 0。
// 递推 a_n = Σ_{i=1..k} c_i a_{n-i} 时:Q(x) = 1 - Σ c_i x^i,P(x) = Σ_{i<k} a_i x^i。
// 复杂度 O(k log k log n)。依赖:mul(poly, poly)(本目录 ntt.cpp 的卷积)。

using poly = vector<mint>;
mint bostan_mori(poly P, poly Q, ll n) {
    // 平凡情形:Q 是常数(deg Q = 0)时 P/Q 只有第 0 项非零。
    // 另外 mul() 对长度 1 的输入会返回空向量,所以这里必须先短路。
    if((int) Q.size() <= 1) return n == 0 ? (P.empty() ? mint(0) : P[0]) / Q[0] : mint(0);
    // 高于 deg Q - 1 的 P 项对第 n 项没有影响,直接截掉(不截会让卷积越算越大)
    if((int) P.size() > (int) Q.size() - 1) P.resize(Q.size() - 1);
    while (n) {
        poly Qm = Q;
        For(i, 1, (int) Qm.size() - 1) if(i & 1) Qm[i] = -Qm[i]; // Q(-x)
        poly U = mul(P, Qm), V = mul(Q, Qm);
        P.clear();
        for(int i = n & 1; i < (int) U.size(); i += 2) P.push_back(U[i]); // 取偶数/奇数项
        Q.clear();
        for(int i = 0; i < (int) V.size(); i += 2) Q.push_back(V[i]);
        n >>= 1;
    }
    return P[0] / Q[0];
}
