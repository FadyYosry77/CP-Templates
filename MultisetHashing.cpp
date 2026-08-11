#include <bits/stdc++.h>
using namespace std;

using ll = long long;
using ull = unsigned long long;

struct MultiHash {
    struct H {
        ull a = 0, b = 0;
        H operator+(H o) const { return {a + o.a, b + o.b}; }
        H operator-(H o) const { return {a - o.a, b - o.b}; }

        H &operator+=(H o) {
            a += o.a, b += o.b;
            return *this;
        }

        H &operator-=(H o) {
            a -= o.a, b -= o.b;
            return *this;
        }

        bool operator==(H o) const { return a == o.a && b == o.b; }
    };

    H h;
    ll sz = 0;
    vector<H> pref, val;

    static ull mix(ull x) {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    static H hash(ll x) {
        static ull s = chrono::steady_clock::now().time_since_epoch().count();
        return {mix(x + s), mix(x + s + 0x123456789abcdefULL)};
    }

    // Add one value to a multiset. O(1)
    void add(ll x) { h += hash(x), sz++; }

    // Remove one value from a multiset. O(1)
    void erase(ll x) { h -= hash(x), sz--; }

    // Clear current multiset. O(1)
    void clear() { h = {}, sz = 0; }

    // Merge two multisets. O(1)
    void add(MultiHash&o) { h += o.h, sz += o.sz; }

    // Remove another multiset. O(1)
    void erase(MultiHash&o) { h -= o.h, sz -= o.sz; }

    // Check if two multisets are equal. O(1)
    bool operator==(MultiHash&o) { return sz == o.sz && h == o.h; }

    // Build prefix multiset hashes. O(n)
    void build(vector<ll>&a) {
        pref.assign(a.size() + 1, {});
        for(ll i = 0;i < a.size();i++) pref[i + 1] = pref[i] + hash(a[i]);
    }

    // Hash of subarray [l,r]. O(1)
    H get(ll l, ll r) { return pref[r + 1] - pref[l]; }

    // Check if two ranges are permutations of each other. O(1)
    bool same(ll l1, ll r1, ll l2, ll r2) {
        return r1 - l1 == r2 - l2 && get(l1, r1) == get(l2, r2);
    }

    // Hash an arbitrary vector as a multiset. O(n)
    static H get(vector<ll>&a) {
        H h;
        for(auto x:a) h += hash(x);
        return h;
    }

    // Precompute hashes of {1..i}. O(mx)
    void build_values(ll mx) {
        val.assign(mx + 1, {});
        for(ll i = 1;i <= mx;i++) val[i] = val[i - 1] + hash(i);
    }

    // Hash of multiset {l,l+1,...,r}. O(1)
    H values(ll l, ll r) { return val[r] - val[l - 1]; }

    // Check if a[l..r] is a permutation of {1..len}. O(1)
    bool perm(ll l, ll r) {
        ll len = r - l + 1;
        return len < val.size() && get(l, r) == val[len];
    }

    // Check if a[l..r] is a permutation of {L..R}. O(1)
    bool perm(ll l, ll r, ll L, ll R) {
        return r - l == R - L && get(l, r) == values(L, R);
    }

    // Check if a[l..r] is a permutation of a fixed pattern. O(1)
    bool perm(ll l, ll r, H target, ll len) {
        return r - l + 1 == len && get(l, r) == target;
    }

    // Check if two whole arrays are permutations. O(n)
    static bool same(vector<ll>&a, vector<ll>&b) {
        return a.size() == b.size() && get(a) == get(b);
    }

    // Check whether replacing x by y makes two multisets equal. O(1)
    static bool replace(H a, H b, ll x, ll y) {
        return a - hash(x) + hash(y) == b;
    }

    // Count windows that are permutations of pattern. O(n+m)
    static ll count(vector<ll>&a, vector<ll>&p) {
        ll n = a.size(), m = p.size(), ans = 0;
        if(m > n) return 0;

        H cur, target = get(p);
        for(ll i = 0;i < m;i++) cur += hash(a[i]);

        ans += cur == target;

        for(ll i = m;i < n;i++) {
            cur += hash(a[i]), cur -= hash(a[i - m]);
            ans += cur == target;
        }
        return ans;
    }

    // Find starts of all windows that are permutations of pattern. O(n+m)
    static vector<ll> find(vector<ll>&a, vector<ll>&p) {
        ll n = a.size(), m = p.size();
        vector<ll> ans;
        if(m > n) return ans;

        H cur, target = get(p);
        for(ll i = 0;i < m;i++) cur += hash(a[i]);

        if(cur == target) ans.push_back(0);

        for(ll i = m;i < n;i++) {
            cur += hash(a[i]), cur -= hash(a[i - m]);
            if(cur == target) ans.push_back(i - m + 1);
        }
        return ans;
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Dynamic multiset
    MultiHash a, b;
    a.add(1), a.add(2), a.add(2), a.add(5);
    b.add(5), b.add(2), b.add(1), b.add(2);
    cout << (a == b) << '\n';

    // Two ranges are permutations
    vector<ll> v = {1, 2, 3, 4, 4, 2, 1, 3};
    MultiHash h;
    h.build(v);
    cout << h.same(0, 3, 4, 7) << '\n';

    // Range is permutation of 1..k
    vector<ll> p = {3, 1, 4, 2, 1, 2, 2, 4};
    h.build(p), h.build_values(100);
    cout << h.perm(0, 3) << '\n';
    cout << h.perm(4, 7) << '\n';

    // Range is permutation of L..R
    vector<ll> q = {7, 5, 6, 12, 10, 11};
    h.build(q), h.build_values(100);
    cout << h.perm(0, 2, 5, 7) << '\n';
    cout << h.perm(3, 5, 10, 12) << '\n';

    // Fixed target permutation
    vector<ll> x = {9, 4, 7, 7, 9, 4}, target = {4, 7, 9};
    h.build(x);
    auto th = MultiHash::get(target);
    cout << h.perm(0, 2, th, target.size()) << '\n';

    // Whole arrays are permutations
    vector<ll> x1 = {1, 2, 2, 5}, x2 = {2, 5, 1, 2};
    cout << MultiHash::same(x1, x2) << '\n';

    // Count permutation windows
    vector<ll> text = {1, 2, 3, 2, 1, 3, 4}, pat = {1, 2, 3};
    cout << MultiHash::count(text, pat) << '\n';

    // Find permutation windows
    for(auto i:MultiHash::find(text, pat)) cout << i << ' ';
    cout << '\n';
}
