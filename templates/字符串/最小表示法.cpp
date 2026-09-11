// 字符串 · 循环同构最小表示法

// a : vector<auto>
// 返回结果：p 表示 a[p..n-1] 最小。
int minrep(auto &a) {
    int n = a.size(), i = 0, j = 1, k = 0;
    while(i < n && j < n && k < n) {
        if(a[(i + k) % n] == a[(j + k) % n]) {
            ++k;
        } else {
            if(a[(i + k) % n] > a[(j + k) % n]) swap(i, j);
            j += k + 1, j += i == j, k = 0;
        }
    }
    return i;
}
// rotate(a.begin(), a.begin() + minrep(a), a.end());