
#include <bits/stdc++.h>

using namespace std;
typedef vector<int> vi;
typedef vector<vi> vv;
typedef long long ll;
typedef unsigned long long ull;
const int N = 4.4e6 + 15, mod = 998244353;
int F[N], X[N], Y[N], W[N];
ull Z[N];
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
	int i, j, w;
	for(i = 0; i <= n; i++)
	{
		for(j = 1; j < (1 << i); j++)
			rev[(1 << i) + j] = rev[(1 << (i - 1)) + (j >> 1)] | ((j & 1) << (i - 1));
		w = qpow(3, (mod - 1) >> i);
		omg[1 << i] = 1;
		for(j = 1; j < (1 << i); j++)
			omg[(1 << i) + j] = (ll)omg[(1 << i) + j - 1] * w % mod;
	}
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
		memset(X, 0, 4 << l);
		memset(Y, 0, 4 << l);
		for(i = 0; i < f.c.size(); i++)
			X[i] = f.c[i];
		for(i = 0; i < g.c.size(); i++)
			Y[i] = g.c[i];
		ntt(X, l, 0);
		ntt(Y, l, 0);
		for(i = 0; i < (1 << l); i++)
			X[i] = (ll)X[i] * Y[i] % mod;
		ntt(X, l, 1);
		h.c.resize(f.c.size() + g.c.size() - 1);
		for(i = 0; i < h.c.size(); i++)
			h.c[i] = X[i];
		return h;
	}
};
vv solve(vi F, vv G, int n, int m) // [x^[0, n)y^(-m, 0]], [x^0]G = 1
{
	int i, j, S, T, a, b, d = 2 * n - 1;
	vv ans;
	if(n == 1)
	{
		F.resize(m);
		ans.resize(m);
		for(i = 0; i < m; i++)
			ans[i] = {F[i]};
		return ans;
	}
	poly A, B, C;
	vv H;
	S = m * d + n;
	A.c.resize(S);
	B.c.resize(S);
	for(i = 0; i <= m; i++)
		for(j = 0; j < n; j++)
		{
			A.c[i * d + j] = G[i][j];
			B.c[i * d + j] = (j & 1 ? mod - G[i][j] : G[i][j]);
		}
	A = A * B;
	a = (n + 1) / 2;
	b = 2 * m;
	H.resize(b + 1);
	for(i = 0; i <= b; i++)
	{
		H[i].resize(a);
		for(j = 0; j < a; j++)
			H[i][j] = A.c[i * d + 2 * j];
	}
	auto P = solve(F, H, a, b);
	T = (b - 1) * d + n;
	C.c.resize(T);
	for(i = 0; i < b; i++)
		for(j = 0; j < a; j++)
			C.c[i * d + 2 * j] = P[b - 1 - i][j];
	C = C * B;
	ans.resize(m);
	for(i = 0; i < m; i++)
	{
		ans[i].resize(n);
		for(j = 0; j < n; j++)
			ans[i][j] = C.c[(b - 1 - i) * d + j];
	}
	return ans;
}
poly comp(poly F, poly G, int n) // [x^0]G = 0
{
	int i;
	vv Gp;
	Gp.resize(2);
	Gp[0].resize(n);
	Gp[1].resize(n);
	Gp[0][0] = 1;
	for(i = 1; i < n; i++)
		Gp[1][i] = mod - G.c[i];
	auto P = solve(F.c, Gp, n, 1);
	return poly(P[0]);
}
int main()
{
	prep(22);
	int n, i;
	poly F, G;
	scanf("%d", &n);
	F.c.resize(n);
	G.c.resize(n);
	for(i = 0; i < n; i++)
		scanf("%d", &F.c[i]);
	for(i = 0; i < n; i++)
		scanf("%d", &G.c[i]);
	auto H = comp(F, G, n);
	for(i = 0; i < n; i++)
		printf("%d ", H.c[i]);
	printf("\n");
	return 0;
}
