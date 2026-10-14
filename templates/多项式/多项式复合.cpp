// comp(): 多项式 · 多项式复合 (F(G(x)) mod x^n, O(n log n))

// 系数低次在前且在 [0, mod)；使用前须 prep(22)，覆盖内部 NTT 长度。
// 两份复合模板独立使用；静态缓冲约 168 MiB，不含临时向量。

typedef vector<int> vi;
typedef vector<vi> vv;
typedef unsigned long long ull;
const int N = 4.4e6 + 15, mod = 998244353;
int X[N], Y[N];
int rev[N * 2], omg[N * 2];
// 预处理 NTT 的 rev/omg 表；要求 n >= 0（按 2^n 长度上界）；用 NTT 前必须先调用一次。
void prep(int n)
{
	int i, j, w;
	for(i = 0; i <= n; i++)
	{
		for(j = 1; j < (1 << i); j++)
			rev[(1 << i) + j] = rev[(1 << (i - 1)) + (j >> 1)] | ((j & 1) << (i - 1));
		w = ksm(3, (mod - 1) >> i, mod);
		omg[1 << i] = 1;
		for(j = 1; j < (1 << i); j++)
			omg[(1 << i) + j] = (ll)omg[(1 << i) + j - 1] * w % mod;
	}
}
ull tmp[N * 2];
// 原地做长度 2^n 的 NTT；要求 A 有 >= 2^n 个元素、已调用 prep(>= n)；t=false 正变换、t=true 逆变换。
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
struct poly
{
	vector<int> c;
	// 默认构造：空多项式。
	poly() = default;
	// 由系数向量构造；d 为低次在前的系数。
	poly(vector<int> d) : c(move(d)) {}
	// 模 mod 卷积；任一为空返回空；返回长度 |f|+|g|-1 的乘积。
	friend poly operator * (poly f, poly g)
	{
		if(f.c.empty() || g.c.empty())
			return poly();
		poly h;
		int i, l = 0;
		while((f.c.size() + g.c.size()) >> l)
			l++;
		memset(X, 0, (1 << l) * sizeof(int));
		memset(Y, 0, (1 << l) * sizeof(int));
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
// 复合的递归分治核心，提取双变量级数的 [x^[0,n)) y^(-m,0] 系数；要求 G 的 x 常数项为 1；返回分治合并结果。
vv solve(vi F, vv G, int n, int m)
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
// 求 F(G(x)) mod x^n；要求 n>=1、G[0]=0（实现忽略 G[0]）。
// F、G 按值传入，低阶不足按零补齐；返回恰好 n 项，不改调用方输入。
// 须先 prep(k) 使预处理覆盖内部 NTT 长度。
poly comp(poly F, poly G, int n)
{
	int i;
	vv Gp(2, vi(n)); // 2×n:第一行 1,第二行 -G
	Gp[0][0] = 1;
	if((int) G.c.size() < n) G.c.resize(n);
	if((int) F.c.size() < n) F.c.resize(n);
	for(i = 1; i < n; i++)
		Gp[1][i] = mod - G.c[i];
	auto P = solve(F.c, Gp, n, 1);
	return poly(P[0]);
}
