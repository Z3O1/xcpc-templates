// CRT(): 数论 · 中国剩余定理 (合并两个同余方程)

// 把 x ≡ a1 (mod p1)、x ≡ a2 (mod p2) 合并成 x ≡ CRT (mod lcm(p1, p2)),
// 返回 [0, lcm) 内的那个解;无解返回 -1。要求 p1, p2 > 0。
ll CRT(ll a1, ll p1, ll a2, ll p2) {
    ll x, y, g = exgcd(p1, p2, x, y);
    if((a2 - a1) % g) return -1;
    ll m1 = p1 / g, m2 = p2 / g;
    // p1 * t ≡ a2 - a1 (mod p2) 约去 g 后:p1 / g * t ≡ (a2 - a1) / g (mod m2)
    ll k = (a2 - a1) / g % m2 * (x % m2 + m2) % m2;
    ll ans = ((i128) k * p1 + a1) % (m1 * p2); // mod lcm 得到 [0, lcm) 内的解
    return (ans + m1 * p2) % (m1 * p2);
}
// 多方程:ans = 0, mod = 1;
//   For(i, 1, n) {
//       ans = CRT(ans, mod, a[i], p[i]);
//       if(ans < 0) 无解;
//       mod = mod / gcd(mod, p[i]) * p[i];
//   }
