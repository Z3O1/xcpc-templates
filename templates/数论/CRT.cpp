// CRT(): 数论 · 中国剩余定理 (合并两个同余方程)

// 合并两个同余方程；返回 [0,lcm(p1,p2)) 内的解，无解返回 -1。
// 要求 0<=a1<p1、0<=a2<p2，且模乘及加法的所有中间值均在 ll 内。
// 统一安全条件为 1<=p1,p2<=2^31；更大模数需另行核验溢出。
ll CRT(ll a1, ll p1, ll a2, ll p2) {
    ll x, y, g = exgcd(p1, p2, x, y);
    if((a2 - a1) % g) return -1;
    ll m1 = p1 / g, m2 = p2 / g, L = m1 * p2; // L = lcm(p1, p2)
    ll k = (a2 - a1) / g % m2 * (x % m2 + m2) % m2;
    ll ans = ((i128) k * p1 + a1) % L; // mod lcm 得到 [0, lcm) 内的解
    return (ans + L) % L;
}
// 多方程:ans = 0, mod = 1;
//   For(i, 1, n) {
//       ans = CRT(ans, mod, a[i], p[i]);
//       if(ans < 0) 无解;
//       mod = mod / gcd(mod, p[i]) * p[i];
//   }