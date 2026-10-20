// minrep(): 字符串 · 循环同构最小表示法

// O(n) 求最小循环表示的起点（0-based，空序列返回 0）。
template <class T> int minrep(const T &a) {
    int n = a.size(), i = 0, j = 1, k = 0;
    while(i < n && j < n && k < n) {
        if(a[(i + k) % n] == a[(j + k) % n]) ++k;
        else {
            if(a[(i + k) % n] > a[(j + k) % n]) swap(i, j);
            j += k + 1, j += (j == i), k = 0;
        }
    }
    return i;
}
