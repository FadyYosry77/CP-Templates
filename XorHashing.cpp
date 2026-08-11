#include <bits/stdc++.h>
using namespace std;

using ll = long long;
using ull = unsigned long long;

struct XorHash {
    struct H {
        ull a = 0, b = 0;
        H operator^(H o) const { return {a ^ o.a, b ^ o.b}; }

        H &operator^=(H o) {
            a ^= o.a, b ^= o.b;
            return *this;
        }

        bool operator==(H o) const { return a == o.a && b == o.b; }
    };

    H h;
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

    // Toggle one value. Two equal occurrences cancel. O(1)
    void toggle(ll x) { h ^= hash(x); }

    // Clear current XOR hash. O(1)
    void clear() { h = {}; }

    // Merge parity information of two collections. O(1)
    void add(XorHash&o) { h ^= o.h; }

    // Current XOR hash. O(1)
    H get() { return h; }

    // Build prefix XOR hashes. O(n)
    void build(vector<ll>&a) {
        pref.assign(a.size() + 1, {});
        for(ll i = 0;i < a.size();i++) pref[i + 1] = pref[i] ^ hash(a[i]);
    }

    // XOR hash of a[l..r]. O(1)
    H get(ll l, ll r) {
        return pref[r + 1] ^ pref[l];
    }

    // Check if two ranges have the same odd-frequency elements. O(1)
    bool same(ll l1, ll r1, ll l2, ll r2) {
        return get(l1, r1) == get(l2, r2);
    }

    // Check if every value occurs an even number of times. O(1)
    bool all_even(ll l, ll r) {
        return get(l, r) == H{};
    }

    // Hash a vector considering only frequency parity. O(n)
    static H get(vector<ll>&a) {
        H h;
        for(auto x:a) h ^= hash(x);
        return h;
    }

    // Check if two arrays have the same frequency parity. O(n)
    static bool same(vector<ll>&a, vector<ll>&b) {
        return get(a) == get(b);
    }

    // Precompute XOR of {1,2,...,i}. O(mx)
    void build_values(ll mx) {
        val.assign(mx + 1, {});
        for(ll i = 1;i <= mx;i++) val[i] = val[i - 1] ^ hash(i);
    }

    // XOR hash of {L,L+1,...,R}. O(1)
    H values(ll L, ll R) {
        return val[R] ^ val[L - 1];
    }

    // Check parity signature against {1,...,len}. O(1)
    bool perm(ll l, ll r) {
        ll len = r - l + 1;
        return len < val.size() && get(l, r) == val[len];
    }

    // Check parity signature against {L,...,R}. O(1)
    bool perm(ll l, ll r, ll L, ll R) {
        return r - l == R - L && get(l, r) == values(L, R);
    }

    // Check if range parity equals a fixed target. O(1)
    bool same(ll l, ll r, H target) {
        return get(l, r) == target;
    }

    // Toggle x hypothetically without modifying original hash. O(1)
    static H with_toggle(H h, ll x) {
        return h ^ hash(x);
    }

    // Replace x with y in parity hash. O(1)
    static H replace(H h, ll x, ll y) {
        return h ^ hash(x) ^ hash(y);
    }

    // Check if replacing x by y makes two parity hashes equal. O(1)
    static bool replace(H a, H b, ll x, ll y) {
        return (a ^ hash(x) ^ hash(y)) == b;
    }

    // Count windows having same frequency parity as pattern. O(n+m)
    static ll count(vector<ll>&a, vector<ll>&p) {
        ll n = a.size(), m = p.size(), ans = 0;
        if(m > n) return 0;

        H cur, target = get(p);
        for(ll i = 0;i < m;i++) cur ^= hash(a[i]);

        ans += cur == target;

        for(ll i = m;i < n;i++) {
            cur ^= hash(a[i]), cur ^= hash(a[i - m]);
            ans += cur == target;
        }
        return ans;
    }

    // Find windows having same frequency parity as pattern. O(n+m)
    static vector<ll> find(vector<ll>&a, vector<ll>&p) {
        ll n = a.size(), m = p.size();
        vector<ll> ans;
        if(m > n) return ans;

        H cur, target = get(p);
        for(ll i = 0;i < m;i++) cur ^= hash(a[i]);

        if(cur == target) ans.push_back(0);

        for(ll i = m;i < n;i++) {
            cur ^= hash(a[i]), cur ^= hash(a[i - m]);
            if(cur == target) ans.push_back(i - m + 1);
        }
        return ans;
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Dynamic toggling
    XorHash a;
    a.toggle(1), a.toggle(2), a.toggle(1);
    cout << (a.get() == XorHash::hash(2)) << '\n';

    // Same odd-frequency elements
    vector<ll> v = {1, 2, 2, 3, 1, 3};
    XorHash h;
    h.build(v);

    cout << h.same(0, 2, 3, 5) << '\n';

    // Every frequency even
    vector<ll> x = {1, 2, 1, 2, 3, 3};
    h.build(x);

    cout << h.all_even(0, 5) << '\n';

    // Prefix parity queries
    vector<ll> q = {1, 2, 3, 2, 1};
    h.build(q);

    cout << (h.get(0, 4) == XorHash::hash(3)) << '\n';

    // Same parity between whole arrays
    vector<ll> x1 = {1, 1, 2, 3};
    vector<ll> x2 = {2, 3};

    cout << XorHash::same(x1, x2) << '\n';

    // Sliding windows with same parity as pattern
    vector<ll> text = {1, 2, 3, 2, 1, 3, 4};
    vector<ll> pat = {1, 2, 3};

    cout << XorHash::count(text, pat) << '\n';

    for(auto i:XorHash::find(text, pat)) cout << i << ' ';
    cout << '\n';
}
