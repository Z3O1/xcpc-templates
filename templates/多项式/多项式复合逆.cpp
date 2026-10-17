// comp_inv(): 多项式 · 多项式复合逆 (复合逆, O(n log n))

// 系数低次在前且在 [0,mod)；须先 prep(k)，0<=k<=22 且覆盖所有内部 NTT 长度。
// 两份复合模板独立使用；静态缓冲约 202 MiB，不含临时向量。

typedef vector<int> vi;
typedef vector<vi> vv;
typedef unsigned long long ull;
const int N = 4.4e6 + 15, mod = 998244353;
int X[N], fac[N], inv[N], Inv[N];
int rev[N * 2], omg[N * 2];
// 预处理长度不超过 2^n 的 NTT 及阶乘/逆元表（0<=n<=22）。
void prep(int n)
{
	int i, j, w, m = N - 5;
	for(i = 0; i <= n; i++)
	{
		for(j = 1; j < (1 << i); j++)
			rev[(1 << i) + j] = rev[(1 << (i - 1)) + (j >> 1)] | ((j & 1) << (i - 1));
		w = ksm(3, (mod - 1) >> i, mod);
		omg[1 << i] = 1;
		for(j = 1; j < (1 << i); j++)
			omg[(1 << i) + j] = (ll)omg[(1 << i) + j - 1] * w % mod;
	}
	fac[0] = 1;
	for(i = 1; i <= m; i++)
		fac[i] = (ll)fac[i - 1] * i % mod;
	inv[m] = ksm(fac[m], mod - 2, mod);
	for(i = m; i >= 1; i--)
		inv[i - 1] = (ll)inv[i] * i % mod;
	for(i = 1; i <= m; i++)
		Inv[i] = (ll)inv[i] * fac[i - 1] % mod;
}
ull tmp[N * 2];
// 原地变换 A[0..2^n)；须已 prep(>=n)，t=false 正变换、t=true 逆变换。
inline void ntt(int *A, int n, bool t)
{
	int i, j, k, q, x;
	for(i = 0; i < (1 << n); i++)
		tmp[i] = A[i];
	for(i = 0; i < (1 << n); i++)
		if(rev[(1 << n) + i] > i)
			swap(tmp[rev[(1 << n) + i]], tmp[i]);
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < (1 << n); j += (2 << i))
			for(k = 0; k < (1 << i); k++)
			{
				x = tmp[j + k + (1 << i)] * omg[(2 << i) + k] % mod;
				tmp[j + k + (1 << i)] = tmp[j + k] + mod - x;
				tmp[j + k] = tmp[j + k] + x;
			}
		if(i == 17 || i == n - 1)
			for(j = 0; j < (1 << n); j++)
				tmp[j] %= mod;
	}
	for(i = 0; i < (1 << n); i++)
		A[i] = tmp[i];
	if(t)
	{
		for(i = 1; i < (1 << n); i++)
			if(i < (1 << n) - i)
				swap(A[i], A[(1 << n) - i]);
		q = ksm(1 << n, mod - 2, mod);
		for(i = 0; i < (1 << n); i++)
			A[i] = (ll)A[i] * q % mod;
	}
}
// 将 F 补到 2^n 项后变换；|F|<=2^n，其余同 ntt，返回新向量。
vi pv(vi F, int n, bool t)
{
	int i;
	F.resize(1 << n);
	for(i = 0; i < (1 << n); i++)
		X[i] = F[i];
	ntt(X, n, t);
	for(i = 0; i < (1 << n); i++)
		F[i] = X[i];
	return F;
}
struct poly
{
	vector<int> c;
	poly() = default;
	poly(vector<int> d) : c(move(d)) {}
	// 卷积保留 |f|+|g|-1 项，任一为空则返回空。
	friend poly operator * (poly f, poly g)
	{
		if(f.c.empty() || g.c.empty())
			return poly();
		poly h;
		int i, l = 0, len = f.c.size() + g.c.size() - 1;
		while((f.c.size() + g.c.size()) >> l)
			l++;
		f.c = pv(f.c, l, 0);
		g.c = pv(g.c, l, 0);
		for(i = 0; i < (1 << l); i++)
			X[i] = (ll)f.c[i] * g.c[i] % mod;
		ntt(X, l, 1);
		h.c.resize(len);
		for(i = 0; i < h.c.size(); i++)
			h.c[i] = X[i];
		return h;
	}
};
// 返回 [x^(n-1)]1/(1-yF(x)) 的前 n 项 y 系数；n>=1、|F.c|>=n，仅供复合逆使用。
poly pw_pj(poly F, int n)
{
	int i, j, S, R, m = 1, d = n - 1, l, r;
	vv A, B;
	vi pA, pB, pC;
	A.resize(2);
	A[0].resize(n);
	A[0][0] = 1;
	A[1].resize(n);
	B = A;
	for(i = 0; i < n; i++)
		B[1][i] = mod - F.c[i];
	while(d)
	{
		R = 2 * d + 1;
		S = m * R + d + 1;
		l = 0;
		while((S - 1) >> l)
			l++;
		l++;
		pA = pB = pC = {};
		pA.resize(S);
		pB.resize(S);
		pC.resize(S);
		for(i = 0; i <= m; i++)
			for(j = 0; j <= d; j++)
			{
				pA[i * R + j] = A[i][j];
				pB[i * R + j] = B[i][j];
				pC[i * R + j] = (j & 1 ? mod - B[i][j] : B[i][j]);
			}
		pA = pv(pA, l, 0);
		pB = pv(pB, l, 0);
		pC = pv(pC, l, 0);
		for(i = 0; i < (1 << l); i++)
		{
			pA[i] = (ll)pA[i] * pC[i] % mod;
			pB[i] = (ll)pB[i] * pC[i] % mod;
		}
		pA = pv(pA, l, 1);
		pB = pv(pB, l, 1);
		r = d & 1;
		d /= 2;
		m *= 2;
		A.resize(m + 1);
		B.resize(m + 1);
		for(i = 0; i <= m; i++)
		{
			A[i].resize(d + 1);
			B[i].resize(d + 1);
			for(j = 0; j <= d; j++)
			{
				A[i][j] = pA[i * R + j * 2 + r];
				B[i][j] = pB[i * R + j * 2];
			}
		}
	}
	F.c.resize(n); // 截去高次项，避免随后 reverse 把它们倒到低次位置。
	for(i = 0; i < n; i++)
		F.c[i] = A[i][0];
	return F;
}
// 返回 n 项的复合逆 G，使 F(G(x)) ≡ x mod x^n；n>=1、F[0]=0、F[1]!=0。
// 不足项补零，不改输入；即使 n=1 也须提供非零一次项。
poly comp_inv(poly F, int n)
{
	if((int) F.c.size() < std::max(n, 2)) F.c.resize(std::max(n, 2));
	int i, r, p = ksm(F.c[1], mod - 2, mod), c;
	for(i = 0, c = 1; i < n; i++, c = (ll)c * p % mod)
		F.c[i] = (ll)F.c[i] * c % mod;
	poly G;
	F = pw_pj(F, n);
	reverse(F.c.begin(), F.c.end());
	for(i = 0; i < n; i++)
		F.c[i] = (ll)F.c[i] * Inv[n - 1 - i] % mod * (n - 1) % mod;
	r = mod - Inv[n - 1];
	G.c = {1};
	while(r)
	{
		if(r & 1)
		{
			G = G * F;
			G.c.resize(n);
		}
		F = F * F;
		F.c.resize(n);
		r >>= 1;
	}
	G.c.resize(n);
	for(i = n - 1; i > 0; i--)
		G.c[i] = (ll)G.c[i - 1] * p % mod;
	G.c[0] = 0;
	return G;
}
