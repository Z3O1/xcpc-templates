
#include <bits/stdc++.h>

using namespace std;
typedef vector<int> vi;
typedef vector<vi> vv;
typedef long long ll;
typedef unsigned long long ull;
const int N = 4.4e6 + 15, mod = 998244353;
int F[N], X[N], fac[N], inv[N], Inv[N];
int rev[N * 2], omg[N * 2];
int qpow(int x, int y)
{
	int ans = 1;
	while(y)
	{
		if(y & 1)
			ans = (ll)ans * x % mod;
		x = (ll)x * x % mod;
		y >>= 1;
	}
	return ans;
}
void prep(int n)
{
	int i, j, w, m = N - 5;
	for(i = 0; i <= n; i++)
	{
		for(j = 1; j < (1 << i); j++)
			rev[(1 << i) + j] = rev[(1 << (i - 1)) + (j >> 1)] | ((j & 1) << (i - 1));
		w = qpow(3, (mod - 1) >> i);
		omg[1 << i] = 1;
		for(j = 1; j < (1 << i); j++)
			omg[(1 << i) + j] = (ll)omg[(1 << i) + j - 1] * w % mod;
	}
	fac[0] = 1;
	for(i = 1; i <= m; i++)
		fac[i] = (ll)fac[i - 1] * i % mod;
	inv[m] = qpow(fac[m], mod - 2);
	for(i = m; i >= 1; i--)
		inv[i - 1] = (ll)inv[i] * i % mod;
	for(i = 1; i <= m; i++)
		Inv[i] = (ll)inv[i] * fac[i - 1] % mod;
}
ull tmp[N * 2];
inline void ntt(int *A, int n, bool t) // 0 DFT 1 IDFT 
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
		q = qpow(1 << n, mod - 2);
		for(i = 0; i < (1 << n); i++)
			A[i] = (ll)A[i] * q % mod;
	}
}
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
	poly()
	{
		c.clear();
	}
	poly(vector<int> d)
	{
		c = d;
	}
	friend poly operator * (poly f, poly g)
	{
		if(f.c.empty() || g.c.empty())
			return poly();
		poly h;
		int i, l = 0;
		while((f.c.size() + g.c.size()) >> l)
			l++;
		f.c = pv(f.c, l, 0);
		g.c = pv(g.c, l, 0);
		for(i = 0; i < (1 << l); i++)
			X[i] = (ll)f.c[i] * g.c[i] % mod;
		ntt(X, l, 1);
		h.c.resize(f.c.size() + g.c.size() - 1);
		for(i = 0; i < h.c.size(); i++)
			h.c[i] = X[i];
		return h;
	}
};
poly pw_pj(poly F, int n) // [x^{n-1}]1/(1-yF)
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
	F.c.resize(n); // 必须截到 n:否则调用方传进来的更长多项式会把高次项留在后面,
	               // 接着 comp_inv 里的 reverse 会把它们倒到低次位置参与运算 → 静默算错
	for(i = 0; i < n; i++)
		F.c[i] = A[i][0];
	return F;
}
// 要求 F.c.size() >= n(n == 1 时至少 2 个系数,内部会读 F.c[1])且 F[1] != 0;
// 返回长度恰为 n 的复合逆。
poly comp_inv(poly F, int n) // [x^n]F^k=k/n[x^{n-k}](G/x)^{-n}
{
	// 补齐系数:下面要读 F.c[1](n == 1 时也得有),pw_pj 也要求 size >= n
	if((int) F.c.size() < std::max(n, 2)) F.c.resize(std::max(n, 2));
	int i, r, p = qpow(F.c[1], mod - 2), c;
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
int main()
{
	prep(22);
	int n, i;
	scanf("%d", &n);
	poly F;
	F.c.resize(n);
	for(i = 0; i < n; i++)
		scanf("%d", &F.c[i]);
	F = comp_inv(F, n);
	for(i = 0; i < n; i++)
		printf("%d ", F.c[i]);
	printf("\n");
	return 0;
}
