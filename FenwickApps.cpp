#include <bits/stdc++.h>
using namespace std;

using ll = long long;
__extension__ typedef __int128 i128;

/*
FENWICK TREE / BINARY INDEXED TREE -- C++17 (also C++20)

QUICK USE
    FenwickTree fw(n);             point add, prefix/range sum
    fw.update(i,v);                a[i] += v (not assignment)
    fw.query(r);                   sum [0,r], and query(-1)==0
    fw.query(l,r);                 sum [l,r]

    RangeAddPointQueryFenwick d(n); d.update(l,r,v), d.query(i)
    RangeFenwickTree rf(n);         rf.update(l,r,v), rf.query(l,r)
    PolynomialFenwickTree pf(n);    add quadratic on a range, range sum
    FenwickMax mx(n);               monotone point chmax, prefix maximum
    FenwickOrderStatisticMultiset   compressed multiset, kth/top-k/sums
    FenwickTree2D                    dense point add + rectangle sum
    SparseFenwickTree2D              predeclared sparse point coordinates

COMMON APPLICATIONS
    inversion_count, count_smaller_after, count_greater_before
    reverse_pairs, count_increasing_triplets, count_crossings
    chord_intersections, josephus_order, remove_by_orders
    lis_length_fenwick, max_weight_increasing_subsequence
    count_increasing_subsequences, lis_length_and_count
    range_count_leq, range_sum_leq, distinct_range_queries
    count_range_sum, count_positive_sum_subarrays
    sliding_window_inversions, sliding_window_median_cost
    nested_interval_counts, subtree/ancestor/HLD helpers
    dominance_2d, offline_rectangle_sums, orthogonal_intersections
    dominance_3d_cdq, earliest_owned_threshold

INTUITION
    Internally index i stores an aggregate for [i-lowbit(i)+1,i], where
    lowbit(i)=i&-i. update uses i+=lowbit(i) to visit larger containing
    intervals; prefix query uses i-=lowbit(i) to partition a prefix.

RECOGNITION
    frequency/rank                 -> compress values, BIT stores counts
    range effects                  -> boundary events / difference BIT
    threshold queries              -> sort by threshold, activate positions
    subtree                        -> Euler tour turns it into an interval
    DP over keys <= x              -> BIT stores sum/max/custom state by key
    2D relation                    -> sweep one coordinate, BIT on the other
    3D relation                    -> CDQ on one dimension, BIT on another
    earliest monotone time         -> parallel binary search + BIT

WARNINGS
    Range kth needs nonnegative frequencies, so prefix sums are monotone.
    Sum/XOR have inverses; min/max do not. FenwickMax cannot delete normally.
    Equal-key processing order decides strict versus non-strict inequalities.
    Dense 2D BIT costs O(n*m). Polynomial expressions may need __int128.
*/

struct FenwickTree {
    using T = ll;
    vector<T> b;

    FenwickTree(int n = 0) : b(n + 1, 0) {
    }

    FenwickTree(const vector<T>&a) { build(a); }

    // Problem:
    //     Build a Fenwick tree from all array values.
    // Fenwick stores:
    //     Sums of its canonical lowbit intervals.
    // Key observation:
    //     Internal cell i propagates once to j=i+lowbit(i).
    // Extra technique:
    //     Linear Fenwick construction.
    // Returns:
    //     Nothing; replaces this tree with the built array.
    // Complexity:
    //     O(n) time and O(n) memory.
    // Notes:
    //     Public indexes remain zero-based.
    void build(const vector<T>&a) {
        int n = a.size();
        b.assign(n + 1, 0);
        for(int i = 1;i <= n;i++) {
            b[i] += a[i - 1];
            int j = i + (i & -i);
            if(j <= n) b[j] += b[i];
        }
    }

    // Problem:
    //     Add v to one array position.
    // Fenwick stores:
    //     Prefix-summable point contributions.
    // Key observation:
    //     Every containing lowbit interval must receive v.
    // Extra technique:
    //     None.
    // Returns:
    //     Nothing; performs a[i] += v, not assignment.
    // Complexity:
    //     O(log n).
    // Notes:
    //     i==size() intentionally performs no writes, useful for r+1 events.
    void update(int i, T v) {
        for(i++;i < (int) b.size();i += i & -i) b[i] += v;
    }

    // Problem:
    //     Query sum a[0..r].
    // Fenwick stores:
    //     Sums of canonical lowbit intervals.
    // Key observation:
    //     Repeatedly removing lowbit partitions the prefix.
    // Extra technique:
    //     None.
    // Returns:
    //     Prefix sum; query(-1) naturally returns 0.
    // Complexity:
    //     O(log n).
    // Notes:
    //     Requires -1 <= r < size().
    T query(int r) const {
        T ret = 0;
        for(r++;r >= 1;r -= r & -r) ret += b[r];
        return ret;
    }

    // Problem:
    //     Query sum a[l..r].
    // Fenwick stores:
    //     Prefix sums implicitly.
    // Key observation:
    //     Range sum is prefix(r)-prefix(l-1).
    // Extra technique:
    //     Group inverse (subtraction).
    // Returns:
    //     Inclusive range sum, or 0 for l>r.
    // Complexity:
    //     O(log n).
    // Notes:
    //     This formula does not generalize to min/max.
    T query(int l, int r) const {
        return l > r ? 0 : query(r) - query(l - 1);
    }

    T total() const { return query((int) b.size() - 2); }

    // Problem:
    //     Find the first index whose prefix sum is at least target.
    // Fenwick stores:
    //     Nonnegative frequencies or weights.
    // Key observation:
    //     Binary lifting walks over whole Fenwick blocks in O(log n).
    // Extra technique:
    //     Order statistics on monotone prefix sums.
    // Returns:
    //     A zero-based index; size() if target>total. target<=0 returns 0.
    // Complexity:
    //     O(log n), not O(log^2 n).
    // Notes:
    //     Invalid if point values can be negative.
    int lower_bound(T target) const {
        int n = b.size() - 1;
        if(target <= 0) return 0;
        int pos = 0, pw = 1;
        T sum = 0;
        while((pw << 1) <= n) pw <<= 1;
        for(;pw;pw >>= 1) {
            int next = pos + pw;
            if(next <= n && sum + b[next] < target)
                sum += b[next], pos = next;
        }
        return pos;
    }

    int kth(T k) const { return lower_bound(k); }
};

// Problem:
//     Apply an assignment a[i]=value while using an additive Fenwick.
// Fenwick stores:
//     Deltas; the separate array stores current point values.
// Key observation:
//     Assignment equals adding value-a[i].
// Extra technique:
//     Mirror array.
// Returns:
//     Nothing.
// Complexity:
//     O(log n).
// Notes:
//     A normal Fenwick does not remember individual values by itself.
void assign_value(FenwickTree&fw, vector<ll>&a, int i, ll value) {
    assert(0 <= i && i < (int)a.size());
    fw.update(i, value - a[i]);
    a[i] = value;
}

struct RangeAddPointQueryFenwick {
    FenwickTree difference;

    RangeAddPointQueryFenwick(int n = 0) : difference(n) {
    }

    RangeAddPointQueryFenwick(const vector<ll>&a) : difference((int) a.size()) {
        for(int i = 0;i < (int) a.size();i++)
            difference.update(i, a[i] - (i ? a[i - 1] : 0));
    }

    // Problem:
    //     Add v to every position in [l,r].
    // Fenwick stores:
    //     Difference-array boundary events.
    // Key observation:
    //     Add +v at l and -v at r+1.
    // Extra technique:
    //     Difference array.
    // Returns:
    //     Nothing.
    // Complexity:
    //     O(log n).
    // Notes:
    //     r+1==n is safely ignored by FenwickTree::update.
    void update(int l, int r, ll v) {
        if(l > r) return;
        difference.update(l, v);
        difference.update(r + 1, -v);
    }

    // Problem:
    //     Read the current value at position i after range additions.
    // Fenwick stores:
    //     Difference-array events.
    // Key observation:
    //     Prefix sum of differences reconstructs a[i].
    // Extra technique:
    //     Difference array.
    // Returns:
    //     Current a[i].
    // Complexity:
    //     O(log n).
    // Notes:
    //     This wrapper does not provide range sums.
    ll query(int i) const { return difference.query(i); }
};

struct RangeFenwickTree {
    FenwickTree slope, intercept;

    RangeFenwickTree(int size = 0) : slope(size), intercept(size) {
    }

    void add_event(int i, ll v) {
        slope.update(i, v);
        intercept.update(i, -v * (i - 1LL));
    }

    // Problem:
    //     Add v to every position in [l,r].
    // Fenwick stores:
    //     Difference d and weighted difference d*(i-1).
    // Key observation:
    //     prefix(r)=(r+1)*sum(d[0..r])-sum(i*d[i]).
    // Extra technique:
    //     Two Fenwick trees and discrete integration.
    // Returns:
    //     Nothing.
    // Complexity:
    //     O(log n).
    // Notes:
    //     External indexes are zero-based.
    void update(int l, int r, ll v) {
        if(l > r) return;
        add_event(l, v);
        add_event(r + 1, -v);
    }

    // Problem:
    //     Query sum a[0..r] after range additions.
    // Fenwick stores:
    //     Difference slopes and integration intercepts.
    // Key observation:
    //     Each d[i] contributes d[i]*(r-i+1) to this prefix.
    // Extra technique:
    //     Two Fenwick trees.
    // Returns:
    //     Prefix sum.
    // Complexity:
    //     O(log n).
    // Notes:
    //     query(-1)==0.
    ll query(int r) const {
        if(r < 0) return 0;
        return slope.query(r) * r + intercept.query(r);
    }

    ll query(int l, int r) const { return l > r ? 0 : query(r) - query(l - 1); }
};

struct PolynomialFenwickTree {
    using T = ll;
    vector<T> b[4];

    PolynomialFenwickTree(int n = 0) {
        for(int t = 0;t < 4;t++) b[t].assign(n + 1, 0);
    }

    T _query(int t, int r) const {
        T ret = 0;
        for(r++;r >= 1;r -= r & -r) ret += b[t][r];
        return ret;
    }

    // Problem:
    //     Query total value on [0,r] after quadratic range additions.
    // Fenwick stores:
    //     Four coefficients of a cubic prefix-sum polynomial, scaled by 6.
    // Key observation:
    //     Sums of 1,x,x^2 are polynomials of degree at most 3.
    // Extra technique:
    //     Discrete integration / Faulhaber formulas.
    // Returns:
    //     Prefix sum.
    // Complexity:
    //     O(log n).
    // Notes:
    //     Intermediate x^3 terms can overflow ll under large constraints.
    T query(int r) const {
        T ret = 0, x = 1;
        for(int t = 0;t < 4;t++) ret += _query(t, r) * x, x *= r;
        return ret / 6;
    }

    T query(int l, int r) const { return l > r ? 0 : query(r) - query(l - 1); }

    void _update(int t, int i, T v) {
        for(i++;i < (int) b[t].size();i += i & -i) b[t][i] += v;
    }

    // Problem:
    //     Add a+b*x+c*x^2 to every integer position x in [l,r].
    // Fenwick stores:
    //     Boundary events for the scaled cubic prefix polynomial.
    // Key observation:
    //     Integrating a quadratic sequence gives a cubic prefix formula.
    // Extra technique:
    //     Four coefficient Fenwicks.
    // Returns:
    //     Nothing.
    // Complexity:
    //     O(log n).
    // Notes:
    //     Interface and formulas preserve the supplied base template.
    void update(ll l, ll r, T a, T bb, T c) {
        if(l > r) return;
        _update(0, l, 6 * a * (-l + 1));
        _update(1, l, 6 * a);
        _update(0, r + 1, -6 * a * (-l + 1) + 6 * a * (r - l + 1));
        _update(1, r + 1, -6 * a);

        _update(0, l, 3 * bb * (-l * l + l));
        _update(1, l, 3 * bb);
        _update(2, l, 3 * bb);
        _update(0, r + 1, -3 * bb * (-l * l + l) + 3 * bb * (r * (r + 1) - (l - 1) * l));
        _update(1, r + 1, -3 * bb);
        _update(2, r + 1, -3 * bb);

        _update(0, l, -c * (l - 1) * l * (2 * l - 1));
        _update(1, l, c);
        _update(2, l, 3 * c);
        _update(3, l, 2 * c);
        _update(0, r + 1, c * (l - 1) * l * (2 * l - 1)
                          + c * (r * (r + 1) * (2 * r + 1) - (l - 1) * l * (2 * l - 1)));
        _update(1, r + 1, -c);
        _update(2, r + 1, -3 * c);
        _update(3, r + 1, -2 * c);
    }

    // Problem:
    //     Add a+b*t+c*t^2 on [l,r], where t=x-l.
    // Fenwick stores:
    //     The equivalent polynomial in absolute coordinate x.
    // Key observation:
    //     Expand (x-l)^2 before calling update.
    // Extra technique:
    //     Coordinate translation.
    // Returns:
    //     Nothing.
    // Complexity:
    //     O(log n).
    // Notes:
    //     Useful for arithmetic and quadratic progressions starting at l.
    void updateRelative(ll l, ll r, T a, T bb, T c) {
        update(l, r, a - bb * l + c * l * l, bb - 2 * c * l, c);
    }

    void add_constant(int l, int r, T v) { updateRelative(l, r, v, 0, 0); }

    void add_arithmetic(int l, int r, T first, T difference) {
        updateRelative(l, r, first, difference, 0);
    }

    void add_quadratic(int l, int r, T a, T bb, T c) {
        updateRelative(l, r, a, bb, c);
    }
};

struct XorFenwickTree {
    vector<ll> b;

    XorFenwickTree(int n = 0) : b(n + 1) {
    }

    void update(int i, ll v) { for(i++;i < (int) b.size();i += i & -i) b[i] ^= v; }

    ll query(int r) const {
        ll ret = 0;
        for(r++;r >= 1;r -= r & -r) ret ^= b[r];
        return ret;
    }

    ll query(int l, int r) const { return l > r ? 0 : query(r) ^ query(l - 1); }
};

struct FenwickMin {
    ll inf;
    vector<ll> b;

    FenwickMin(int size = 0, ll infinity = (1LL << 62)) : inf(infinity), b(size + 1, infinity) {
    }

    // Only chmin updates are valid; there is no deletion or arbitrary range minimum.
    void update(int i, ll value) {
        for(i++;i < (int) b.size();i += i & -i) b[i] = min(b[i], value);
    }

    ll query(int r) const {
        ll ret = inf;
        for(r++;r >= 1;r -= r & -r) ret = min(ret, b[r]);
        return ret;
    }
};

struct FenwickMax {
    ll negInf;
    vector<ll> b;

    FenwickMax(int size = 0, ll negativeInfinity = -(1LL << 62))
        : negInf(negativeInfinity), b(size + 1, negativeInfinity) {
    }

    // Only chmax updates are valid; query is maximum on [0,r].
    void update(int i, ll value) {
        for(i++;i < (int) b.size();i += i & -i) b[i] = max(b[i], value);
    }

    ll query(int r) const {
        ll ret = negInf;
        for(r++;r >= 1;r -= r & -r) ret = max(ret, b[r]);
        return ret;
    }
};

struct LisState {
    int length = 0;
    ll count = 0;
};

LisState merge_lis_state(LisState a, LisState b, ll mod = 0) {
    if(a.length != b.length) return a.length > b.length ? a : b;
    if(!a.length) return {0, 0};
    a.count = mod ? (a.count + b.count) % mod : a.count + b.count;
    return a;
}

struct FenwickLisState {
    vector<LisState> b;
    ll mod;

    FenwickLisState(int size = 0, ll modulus = 0) : b(size + 1), mod(modulus) {
    }

    void update(int i, LisState value) {
        for(i++;i < (int) b.size();i += i & -i) b[i] = merge_lis_state(b[i], value, mod);
    }

    LisState query(int r) const {
        LisState ret;
        for(r++;r >= 1;r -= r & -r) ret = merge_lis_state(ret, b[r], mod);
        return ret;
    }
};

struct ModuloFenwickTree {
    ll mod;
    vector<ll> b;
    ModuloFenwickTree(int size = 0, ll modulus = 1) : mod(modulus), b(size + 1) { assert(mod >= 1); }

    void update(int i, ll v) {
        v %= mod;
        if(v < 0) v += mod;
        for(i++;i < (int) b.size();i += i & -i) b[i] = (b[i] + v) % mod;
    }

    ll query(int r) const {
        ll ret = 0;
        for(r++;r >= 1;r -= r & -r) ret = (ret + b[r]) % mod;
        return ret;
    }

    ll query(int l, int r) const { return (query(r) - query(l - 1) + mod) % mod; }
    ll total() const { return query((int) b.size() - 2); }
};

struct FenwickTree2D {
    int n = 0, m = 0;
    vector<vector<ll> > b;

    FenwickTree2D(int rows = 0, int columns = 0)
        : n(rows), m(columns), b(rows + 1, vector<ll>(columns + 1)) {
    }

    void update(int x, int y, ll v) {
        for(x++;x <= n;x += x & -x)
            for(int j = y + 1;j <= m;j += j & -j) b[x][j] += v;
    }

    ll query(int x, int y) const {
        ll ret = 0;
        for(x++;x >= 1;x -= x & -x)
            for(int j = y + 1;j >= 1;j -= j & -j) ret += b[x][j];
        return ret;
    }

    ll query(int x1, int y1, int x2, int y2) const {
        if(x1 > x2 || y1 > y2) return 0;
        return query(x2, y2) - query(x1 - 1, y2) - query(x2, y1 - 1) + query(x1 - 1, y1 - 1);
    }
};

struct RangeFenwickTree2D {
    int n = 0, m = 0;
    FenwickTree2D one, xWeight, yWeight, xyWeight;

    RangeFenwickTree2D(int rows = 0, int columns = 0)
        : n(rows), m(columns), one(rows, columns), xWeight(rows, columns),
          yWeight(rows, columns), xyWeight(rows, columns) {
    }

    void add_event(int x, int y, ll value) {
        if(x >= n || y >= m) return;
        one.update(x, y, value);
        xWeight.update(x, y, value * x);
        yWeight.update(x, y, value * y);
        xyWeight.update(x, y, value * x * y);
    }

    // Problem:
    //     Add value to every cell of a closed rectangle.
    // Fenwick stores:
    //     Four 2D difference coefficient arrays.
    // Key observation:
    //     Four corner events describe a rectangle; two-dimensional integration needs four moments.
    // Extra technique:
    //     2D difference array and inclusion-exclusion.
    // Returns:
    //     Nothing.
    // Complexity:
    //     O(log n log m).
    // Notes:
    //     Dense memory is 4*O(n*m); do not use it for dimensions near 2e5.
    void update(int x1, int y1, int x2, int y2, ll value) {
        if(x1 > x2 || y1 > y2) return;
        add_event(x1, y1, value);
        add_event(x1, y2 + 1, -value);
        add_event(x2 + 1, y1, -value);
        add_event(x2 + 1, y2 + 1, value);
    }

    // Problem:
    //     Query the sum of all cells in prefix [0..x]*[0..y].
    // Fenwick stores:
    //     Difference count and x/y/xy weighted moments.
    // Key observation:
    //     Each event contributes value*(x-i+1)*(y-j+1).
    // Extra technique:
    //     Two-dimensional discrete integration.
    // Returns:
    //     Prefix rectangle sum.
    // Complexity:
    //     O(log n log m).
    // Notes:
    //     Negative x or y naturally returns zero through the underlying queries.
    ll query(int x, int y) const {
        if(x < 0 || y < 0) return 0;
        return one.query(x, y) * (x + 1LL) * (y + 1LL)
               - xWeight.query(x, y) * (y + 1LL) - yWeight.query(x, y) * (x + 1LL)
               + xyWeight.query(x, y);
    }

    ll query(int x1, int y1, int x2, int y2) const {
        return query(x2, y2) - query(x1 - 1, y2) - query(x2, y1 - 1) + query(x1 - 1, y1 - 1);
    }
};

struct SparseFenwickTree2D {
    int n = 0;
    vector<vector<int> > ys;
    vector<vector<ll> > b;

    // All coordinates that may ever be updated must be supplied in advance.
    SparseFenwickTree2D(int rows = 0, const vector<pair<int, int> >&points = {}) : n(rows), ys(rows + 1) {
        for(auto [x,y]:points) {
            assert(0 <= x && x < n);
            for(int i = x + 1;i <= n;i += i & -i) ys[i].push_back(y);
        }
        b.resize(n + 1);
        for(int i = 1;i <= n;i++) {
            sort(ys[i].begin(), ys[i].end());
            ys[i].erase(unique(ys[i].begin(), ys[i].end()), ys[i].end());
            b[i].assign(ys[i].size() + 1, 0);
        }
    }

    void update(int x, int y, ll value) {
        for(int i = x + 1;i <= n;i += i & -i) {
            int p = lower_bound(ys[i].begin(), ys[i].end(), y) - ys[i].begin() + 1;
            assert(p <= (int)ys[i].size() && ys[i][p-1] == y);
            for(int j = p;j < (int) b[i].size();j += j & -j) b[i][j] += value;
        }
    }

    ll query(int x, int y) const {
        ll ret = 0;
        for(int i = x + 1;i >= 1;i -= i & -i) {
            int p = upper_bound(ys[i].begin(), ys[i].end(), y) - ys[i].begin();
            for(int j = p;j >= 1;j -= j & -j) ret += b[i][j];
        }
        return ret;
    }

    ll query(int x1, int y1, int x2, int y2) const {
        return query(x2, y2) - query(x1 - 1, y2) - query(x2, y1 - 1) + query(x1 - 1, y1 - 1);
    }
};

struct FenwickOrderStatisticMultiset {
    vector<ll> values, frequency;
    FenwickTree countBit, sumBit;
    ll elements = 0, totalSum = 0;

    FenwickOrderStatisticMultiset(vector<ll> coordinates = {}) {
        sort(coordinates.begin(), coordinates.end());
        coordinates.erase(unique(coordinates.begin(), coordinates.end()), coordinates.end());
        values = move(coordinates);
        frequency.assign(values.size(), 0);
        countBit = FenwickTree(values.size());
        sumBit = FenwickTree(values.size());
    }

    int id(ll x) const { return lower_bound(values.begin(), values.end(), x) - values.begin(); }

    void insert(ll x, ll copies = 1) {
        int p = id(x);
        assert(copies >= 0 && p < (int)values.size() && values[p] == x);
        frequency[p] += copies;
        elements += copies;
        totalSum += copies * x;
        countBit.update(p, copies);
        sumBit.update(p, copies * x);
    }

    void erase(ll x, ll copies = 1) {
        int p = id(x);
        assert(copies >= 0 && p < (int)values.size() && values[p] == x);
        assert(frequency[p] >= copies);
        frequency[p] -= copies;
        elements -= copies;
        totalSum -= copies * x;
        countBit.update(p, -copies);
        sumBit.update(p, -copies * x);
    }

    ll count_less(ll x) const {
        int p = lower_bound(values.begin(), values.end(), x) - values.begin();
        return countBit.query(p - 1);
    }

    ll count_leq(ll x) const {
        int p = upper_bound(values.begin(), values.end(), x) - values.begin();
        return countBit.query(p - 1);
    }

    ll sum_less(ll x) const {
        int p = lower_bound(values.begin(), values.end(), x) - values.begin();
        return sumBit.query(p - 1);
    }

    ll sum_leq(ll x) const {
        int p = upper_bound(values.begin(), values.end(), x) - values.begin();
        return sumBit.query(p - 1);
    }

    ll kth(ll k) const {
        assert(1 <= k && k <= elements);
        return values[countBit.kth(k)];
    }

    ll median() const {
        assert(elements);
        return kth((elements + 1) / 2);
    }

    ll sum_k_smallest(ll k) const {
        assert(0 <= k && k <= elements);
        if(!k) return 0;
        int p = countBit.kth(k);
        ll before = countBit.query(p - 1);
        return sumBit.query(p - 1) + (k - before) * values[p];
    }

    ll sum_k_largest(ll k) const {
        assert(0 <= k && k <= elements);
        return totalSum - sum_k_smallest(elements - k);
    }
};

// Problem:
//     Coordinate-compress arbitrary values while preserving their order relation.
// Fenwick stores:
//     Later routines use the returned ranks as BIT indexes.
// Key observation:
//     Only relative order matters for <,<=,>,>= queries.
// Extra technique:
//     Sort and unique.
// Returns:
//     {rank per element, sorted unique values}.
// Complexity:
//     O(n log n).
// Notes:
//     For arbitrary threshold x, lower_bound counts <x and upper_bound counts <=x.
pair<vector<int>, vector<ll> > compress_values(const vector<ll>&a) {
    vector<ll> values = a;
    sort(values.begin(), values.end());
    values.erase(unique(values.begin(), values.end()), values.end());
    vector<int> rank(a.size());
    for(int i = 0;i < (int) a.size();i++)
        rank[i] = lower_bound(values.begin(), values.end(), a[i]) - values.begin();
    return {rank, values};
}

// Problem:
//     Count inversions: pairs i<j with a[i]>a[j].
// Fenwick stores:
//     Frequencies of compressed values already processed.
// Key observation:
//     Previous greater = processed - previous values <= a[i].
// Extra technique:
//     Coordinate compression.
// Returns:
//     Number of strict inversions.
// Complexity:
//     O(n log n).
// Notes:
//     Duplicates are not inversions; use ll because the answer can be O(n^2).
ll inversion_count(const vector<ll>&a) {
    auto [rank,values] = compress_values(a);
    FenwickTree fw(values.size());
    ll ans = 0;
    for(int i = 0;i < (int) a.size();i++) {
        ans += i - fw.query(rank[i]);
        fw.update(rank[i], 1);
    }
    return ans;
}

// Problem:
//     For every i, count j>i with a[j]<a[i].
// Fenwick stores:
//     Frequencies in the processed suffix.
// Key observation:
//     Scan right-to-left and query ranks strictly below rank[i].
// Extra technique:
//     Coordinate compression and reverse sweep.
// Returns:
//     One strict-smaller count per index.
// Complexity:
//     O(n log n).
// Notes:
//     Query rank[i]-1, not rank[i], to exclude equal values.
vector<ll> count_smaller_after(const vector<ll>&a) {
    auto [rank,values] = compress_values(a);
    FenwickTree fw(values.size());
    vector<ll> ans(a.size());
    for(int i = (int) a.size() - 1;i >= 0;i--)
        ans[i] = fw.query(rank[i] - 1), fw.update(rank[i], 1);
    return ans;
}

// Problem:
//     For every i, count j<i with a[j]>a[i].
// Fenwick stores:
//     Frequencies in the processed prefix.
// Key observation:
//     Previous greater = processed - count <= current.
// Extra technique:
//     Coordinate compression.
// Returns:
//     One strict-greater count per index.
// Complexity:
//     O(n log n).
// Notes:
//     Equal values are excluded.
vector<ll> count_greater_before(const vector<ll>&a) {
    auto [rank,values] = compress_values(a);
    FenwickTree fw(values.size());
    vector<ll> ans(a.size());
    for(int i = 0;i < (int) a.size();i++)
        ans[i] = i - fw.query(rank[i]), fw.update(rank[i], 1);
    return ans;
}

// Problem:
//     Count pairs i<j with a[i]>2*a[j].
// Fenwick stores:
//     Frequencies of previous original values.
// Key observation:
//     For each right endpoint, subtract count <= 2*a[j] from processed count.
// Extra technique:
//     Coordinate compression and __int128-safe threshold comparison.
// Returns:
//     Number of reverse pairs.
// Complexity:
//     O(n log n).
// Notes:
//     upper_bound on __int128 is emulated to avoid overflowing 2*a[j].
ll reverse_pairs(const vector<ll>&a) {
    vector<ll> values = a;
    sort(values.begin(), values.end());
    values.erase(unique(values.begin(), values.end()), values.end());
    FenwickTree fw(values.size());
    ll ans = 0;
    for(int i = 0;i < (int) a.size();i++) {
        i128 threshold = (i128) 2 * a[i];
        int lo = 0, hi = values.size();
        while(lo < hi) {
            int mid = (lo + hi) / 2;
            if((i128) values[mid] <= threshold) lo = mid + 1;
            else hi = mid;
        }
        ans += i - fw.query(lo - 1);
        int p = lower_bound(values.begin(), values.end(), a[i]) - values.begin();
        fw.update(p, 1);
    }
    return ans;
}

// Problem:
//     Count triples i<j<k with a[i]<a[j]<a[k].
// Fenwick stores:
//     Prefix/suffix frequencies of compressed values.
// Key observation:
//     For each middle j multiply leftSmaller[j] by rightGreater[j].
// Extra technique:
//     Two directional sweeps.
// Returns:
//     Number of strictly increasing triples.
// Complexity:
//     O(n log n).
// Notes:
//     Duplicates never satisfy either strict inequality.
ll count_increasing_triplets(const vector<ll>&a) {
    auto [rank,values] = compress_values(a);
    int n = a.size();
    vector<ll> left(n), right(n);
    FenwickTree fw(values.size());
    for(int i = 0;i < n;i++) left[i] = fw.query(rank[i] - 1), fw.update(rank[i], 1);
    fw = FenwickTree(values.size());
    for(int i = n - 1;i >= 0;i--) {
        ll seen = n - 1 - i;
        right[i] = seen - fw.query(rank[i]);
        fw.update(rank[i], 1);
    }
    ll ans = 0;
    for(int i = 0;i < n;i++) ans += left[i] * right[i];
    return ans;
}

// Problem:
//     Count triples i<j<k with a[i]>a[j]>a[k].
// Fenwick stores:
//     Prefix/suffix frequencies of compressed values.
// Key observation:
//     For each middle j multiply leftGreater[j] by rightSmaller[j].
// Extra technique:
//     Two directional sweeps.
// Returns:
//     Number of strictly decreasing triples.
// Complexity:
//     O(n log n).
// Notes:
//     Direct rank queries avoid negating LLONG_MIN values.
ll count_decreasing_triplets(const vector<ll>&a) {
    auto [rank,values] = compress_values(a);
    int n = a.size();
    vector<ll> left(n), right(n);
    FenwickTree fw(values.size());
    for(int i = 0;i < n;i++) left[i] = i - fw.query(rank[i]), fw.update(rank[i], 1);
    fw = FenwickTree(values.size());
    for(int i = n - 1;i >= 0;i--) right[i] = fw.query(rank[i] - 1), fw.update(rank[i], 1);
    ll ans = 0;
    for(int i = 0;i < n;i++) ans += left[i] * right[i];
    return ans;
}

// Problem:
//     Count crossings among edges (leftCoordinate,rightCoordinate).
// Fenwick stores:
//     Frequencies of right endpoints from earlier left groups.
// Key observation:
//     For a1<a2, edges cross when b1>b2.
// Extra technique:
//     Sort by first coordinate and batch equal first coordinates.
// Returns:
//     Crossings with strict inequality in both coordinates.
// Complexity:
//     O(n log n).
// Notes:
//     Query a whole equal-a group before inserting it, so equal endpoints do not count.
ll count_crossings(vector<pair<ll, ll> > edges) {
    vector<ll> ys;
    for(auto [x,y]:edges) ys.push_back(y);
    sort(ys.begin(), ys.end());
    ys.erase(unique(ys.begin(), ys.end()), ys.end());
    sort(edges.begin(), edges.end());
    FenwickTree fw(ys.size());
    ll ans = 0, inserted = 0;
    for(int l = 0, r;l < (int) edges.size();l = r) {
        for(r = l;r < (int) edges.size() && edges[r].first == edges[l].first;r++);
        for(int i = l;i < r;i++) {
            int p = lower_bound(ys.begin(), ys.end(), edges[i].second) - ys.begin();
            ans += inserted - fw.query(p);
        }
        for(int i = l;i < r;i++) {
            int p = lower_bound(ys.begin(), ys.end(), edges[i].second) - ys.begin();
            fw.update(p, 1);
            inserted++;
        }
    }
    return ans;
}

// Problem:
//     Count crossing chords from a circular endpoint label sequence.
// Fenwick stores:
//     Active first occurrences indexed by their positions.
// Key observation:
//     Closing chord [first,i] crosses every active start strictly inside it.
// Extra technique:
//     Event activation/deactivation.
// Returns:
//     Number of unordered crossing chord pairs.
// Complexity:
//     O(n log n).
// Notes:
//     Every label must occur exactly twice.
ll chord_intersections(const vector<int>&labels) {
    FenwickTree fw(labels.size());
    unordered_map<int, int> first;
    ll ans = 0;
    for(int i = 0;i < (int) labels.size();i++) {
        auto it = first.find(labels[i]);
        if(it == first.end()) first[labels[i]] = i, fw.update(i, 1);
        else {
            ans += fw.query(it->second + 1, i - 1);
            fw.update(it->second, -1);
            first.erase(it);
        }
    }
    assert(first.empty());
    return ans;
}

// Problem:
//     Produce Josephus elimination order, skipping k-1 alive people each turn.
// Fenwick stores:
//     One at every alive index.
// Key observation:
//     Binary lifting selects the required alive rank directly.
// Extra technique:
//     Fenwick order statistics.
// Returns:
//     Zero-based elimination order.
// Complexity:
//     O(n log n).
// Notes:
//     k is one-based: k=1 removes the current person.
vector<int> josephus_order(int n, ll k) {
    assert(n >= 0 && k >= 1);
    FenwickTree fw(vector<ll>(n, 1));
    vector<int> ans;
    ll pos = 0;
    for(int alive = n;alive;alive--) {
        pos = (pos + k - 1) % alive;
        int at = fw.kth(pos + 1);
        ans.push_back(at);
        fw.update(at, -1);
    }
    return ans;
}

// Problem:
//     Repeatedly remove the k-th currently alive array element.
// Fenwick stores:
//     Alive flags over original indices.
// Key observation:
//     kth converts an alive rank back to an original position.
// Extra technique:
//     Order-statistic Fenwick.
// Returns:
//     Removed values in query order.
// Complexity:
//     O((n+q) log n).
// Notes:
//     Each order is one-based among currently alive elements.
vector<ll> remove_by_orders(const vector<ll>&a, const vector<int>&orders) {
    FenwickTree fw(vector<ll>(a.size(), 1));
    vector<ll> ans;
    for(int k:orders) {
        assert(1 <= k && k <= fw.total());
        int at = fw.kth(k);
        ans.push_back(a[at]);
        fw.update(at, -1);
    }
    return ans;
}

// Problem:
//     Find the length of a strictly increasing subsequence.
// Fenwick stores:
//     Maximum DP length by compressed ending value.
// Key observation:
//     dp[i]=1+max(dp[j]) over j<i and a[j]<a[i].
// Extra technique:
//     Coordinate compression and prefix-maximum DP.
// Returns:
//     LIS length.
// Complexity:
//     O(n log n).
// Notes:
//     Patience sorting is simpler for length only; this pattern generalizes.
int lis_length_fenwick(const vector<ll>&a) {
    auto [rank,values] = compress_values(a);
    FenwickMax fw(values.size(), 0);
    ll ans = 0;
    for(int p:rank) {
        ll cur = 1 + fw.query(p - 1);
        fw.update(p, cur);
        ans = max(ans, cur);
    }
    return ans;
}

// Problem:
//     Maximum total weight of a strictly increasing value subsequence.
// Fenwick stores:
//     Best weight ending at each compressed value.
// Key observation:
//     Extend the best state with ending value strictly smaller than value[i].
// Extra technique:
//     Weighted LIS DP.
// Returns:
//     Maximum weight, allowing the empty subsequence of weight 0.
// Complexity:
//     O(n log n).
// Notes:
//     For a nonempty-required version with negative weights, adjust the base state.
ll max_weight_increasing_subsequence(const vector<ll>&value, const vector<ll>&weight) {
    assert(value.size() == weight.size());
    auto [rank,values] = compress_values(value);
    FenwickMax fw(values.size(), 0);
    ll ans = 0;
    for(int i = 0;i < (int) value.size();i++) {
        ll cur = fw.query(rank[i] - 1) + weight[i];
        fw.update(rank[i], cur);
        ans = max(ans, cur);
    }
    return ans;
}

// Problem:
//     Count strictly increasing subsequences of exactly length k.
// Fenwick stores:
//     For each length, counts indexed by compressed last value.
// Key observation:
//     A length-l subsequence appends a[i] to every length-(l-1) state ending below it.
// Extra technique:
//     DP with k Fenwick trees.
// Returns:
//     Count modulo mod.
// Complexity:
//     O(n*k*log n) time and O(n*k) BIT memory.
// Notes:
//     Update lengths descending so an element is not reused in the same iteration.
ll count_increasing_subsequences(const vector<ll>&a, int k, ll mod) {
    assert(k >= 1 && mod >= 1);
    auto [rank,values] = compress_values(a);
    vector<ModuloFenwickTree> fw(k + 1, ModuloFenwickTree(values.size(), mod));
    for(int p:rank)
        for(int len = k;len >= 1;len--) {
            ll ways = len == 1 ? 1 : fw[len - 1].query(p - 1) % mod;
            fw[len].update(p, ways);
        }
    return fw[k].total() % mod;
}

// Problem:
//     Count all nonempty strictly increasing subsequences.
// Fenwick stores:
//     Number of subsequences by compressed final value.
// Key observation:
//     ways[i]=1+sum of ways ending with a smaller value.
// Extra technique:
//     Fenwick-sum DP.
// Returns:
//     Count modulo mod.
// Complexity:
//     O(n log n).
// Notes:
//     Query rank-1 for strict; query rank for nondecreasing subsequences.
ll count_all_increasing_subsequences(const vector<ll>&a, ll mod) {
    assert(mod >= 1);
    auto [rank,values] = compress_values(a);
    ModuloFenwickTree fw(values.size(), mod);
    ll ans = 0;
    for(int p:rank) {
        ll ways = (1 + fw.query(p - 1)) % mod;
        fw.update(p, ways);
        ans = (ans + ways) % mod;
    }
    return ans;
}

// Problem:
//     Find LIS length and number of LISs.
// Fenwick stores:
//     Pair (best length,number of ways) for each value prefix.
// Key observation:
//     Merge keeps larger length; ties add counts.
// Extra technique:
//     Custom associative Fenwick state.
// Returns:
//     {LIS length,count}, optionally modulo mod when mod>0.
// Complexity:
//     O(n log n).
// Notes:
//     Equal values query only smaller ranks and therefore cannot extend one another.
pair<int, ll> lis_length_and_count(const vector<ll>&a, ll mod = 0) {
    if(a.empty()) return {0, 0};
    auto [rank,values] = compress_values(a);
    FenwickLisState fw(values.size(), mod);
    for(int p:rank) {
        LisState cur = fw.query(p - 1);
        cur = {cur.length + 1, cur.length ? cur.count : 1};
        fw.update(p, cur);
    }
    LisState ans = fw.query(values.size() - 1);
    return {ans.length, ans.count};
}

struct ThresholdRangeQuery {
    int l = 0, r = -1;
    ll x = 0;
};

// Problem:
//     For every query [l,r,x], count positions i in [l,r] with a[i]<=x.
// Fenwick stores:
//     One at positions whose value has been activated.
// Key observation:
//     Sort elements and queries by value x; sweep x upward.
// Extra technique:
//     Offline sorting by threshold.
// Returns:
//     Answers in original query order.
// Complexity:
//     O((n+q) log(n+q)).
// Notes:
//     Activating <=x before answering makes the inequality non-strict.
vector<ll> range_count_leq(const vector<ll>&a, const vector<ThresholdRangeQuery>&queries) {
    vector<int> positions(a.size()), order(queries.size());
    iota(positions.begin(), positions.end(), 0);
    iota(order.begin(), order.end(), 0);
    sort(positions.begin(), positions.end(), [&](int i, int j) { return a[i] < a[j]; });
    sort(order.begin(), order.end(), [&](int i, int j) { return queries[i].x < queries[j].x; });
    FenwickTree fw(a.size());
    vector<ll> ans(queries.size());
    int at = 0;
    for(int id:order) {
        while(at < (int) a.size() && a[positions[at]] <= queries[id].x)
            fw.update(positions[at++], 1);
        ans[id] = fw.query(queries[id].l, queries[id].r);
    }
    return ans;
}

// Problem:
//     For every query [l,r,x], sum a[i] over i in [l,r] with a[i]<=x.
// Fenwick stores:
//     Activated element values at their original positions.
// Key observation:
//     The same threshold sweep works with weights instead of counts.
// Extra technique:
//     Offline activation by value.
// Returns:
//     Weighted answers in query order.
// Complexity:
//     O((n+q) log(n+q)).
// Notes:
//     Count >x is length-count<=x; sum >x is intervalSum-sum<=x.
vector<ll> range_sum_leq(const vector<ll>&a, const vector<ThresholdRangeQuery>&queries) {
    vector<int> positions(a.size()), order(queries.size());
    iota(positions.begin(), positions.end(), 0);
    iota(order.begin(), order.end(), 0);
    sort(positions.begin(), positions.end(), [&](int i, int j) { return a[i] < a[j]; });
    sort(order.begin(), order.end(), [&](int i, int j) { return queries[i].x < queries[j].x; });
    FenwickTree fw(a.size());
    vector<ll> ans(queries.size());
    int at = 0;
    for(int id:order) {
        while(at < (int) a.size() && a[positions[at]] <= queries[id].x) {
            int p = positions[at++];
            fw.update(p, a[p]);
        }
        ans[id] = fw.query(queries[id].l, queries[id].r);
    }
    return ans;
}

// Problem:
//     Answer number of distinct values in each static subarray [l,r].
// Fenwick stores:
//     One only at the latest processed occurrence of each value.
// Key observation:
//     At sweep endpoint r, latest occurrences inside [l,r] correspond one-to-one to distinct values.
// Extra technique:
//     Offline queries by right endpoint and last occurrence events.
// Returns:
//     Answers in original query order.
// Complexity:
//     O((n+q) log n) expected time with unordered_map.
// Notes:
//     Remove the old latest occurrence before adding the new one.
vector<ll> distinct_range_queries(const vector<ll>&a, const vector<pair<int, int> >&queries) {
    int n = a.size();
    vector<vector<pair<int, int> > > byRight(n);
    vector<ll> ans(queries.size());
    for(int id = 0;id < (int) queries.size();id++) {
        auto [l,r] = queries[id];
        assert(0 <= l && l <= r && r < n);
        byRight[r].push_back({l, id});
    }
    unordered_map<ll, int> last;
    FenwickTree fw(n);
    for(int r = 0;r < n;r++) {
        auto it = last.find(a[r]);
        if(it != last.end()) fw.update(it->second, -1);
        fw.update(r, 1);
        last[a[r]] = r;
        for(auto [l,id]:byRight[r]) ans[id] = fw.query(l, r);
    }
    return ans;
}

// Problem:
//     Answer how many values occur exactly once in every static [l,r].
// Fenwick stores:
//     Difference events over possible left endpoints for the current right endpoint.
// Key observation:
//     If the last occurrences are pp<p<i, value is unique exactly for l in (p,i].
// Extra technique:
//     Offline right-endpoint sweep plus range-add/point-query.
// Returns:
//     Answers in original query order.
// Complexity:
//     O((n+q) log n) expected time.
// Notes:
//     This is a different invariant from ordinary distinct counting.
vector<ll> exactly_once_range_queries(const vector<ll>&a,
                                      const vector<pair<int, int> >&queries) {
    int n = a.size();
    vector<vector<pair<int, int> > > byRight(n);
    vector<ll> ans(queries.size());
    for(int id = 0;id < (int) queries.size();id++) {
        auto [l,r] = queries[id];
        assert(0 <= l && l <= r && r < n);
        byRight[r].push_back({l, id});
    }
    unordered_map<ll, pair<int, int> > previous;
    RangeAddPointQueryFenwick fw(n);
    for(int i = 0;i < n;i++) {
        auto [pp,p] = previous.count(a[i]) ? previous[a[i]] : pair<int, int>{-1, -1};
        if(p != -1) fw.update(pp + 1, p, -1);
        fw.update(p + 1, i, 1);
        previous[a[i]] = {p, i};
        for(auto [l,id]:byRight[i]) ans[id] = fw.query(l);
    }
    return ans;
}

// Problem:
//     Count subarrays with sum in inclusive interval [lower,upper].
// Fenwick stores:
//     Frequencies of previous prefix sums.
// Key observation:
//     lower<=pref[j]-pref[i]<=upper iff pref[j]-upper<=pref[i]<=pref[j]-lower.
// Extra technique:
//     Prefix sums and coordinate compression.
// Returns:
//     Number of qualifying subarrays.
// Complexity:
//     O(n log n).
// Notes:
//     lower>upper returns 0; all threshold searches allow arbitrary values.
ll count_range_sum(const vector<ll>&a, ll lower, ll upper) {
    if(lower > upper) return 0;
    vector<ll> prefix(1, 0);
    for(ll x:a) prefix.push_back(prefix.back() + x);
    vector<ll> values = prefix;
    sort(values.begin(), values.end());
    values.erase(unique(values.begin(), values.end()), values.end());
    FenwickTree fw(values.size());
    ll ans = 0;
    for(ll current:prefix) {
        ll lo = current - upper, hi = current - lower;
        int l = lower_bound(values.begin(), values.end(), lo) - values.begin();
        int r = upper_bound(values.begin(), values.end(), hi) - values.begin() - 1;
        ans += fw.query(l, r);
        fw.update(lower_bound(values.begin(), values.end(), current) - values.begin(), 1);
    }
    return ans;
}

// Problem:
//     Count subarrays whose sum is strictly less than k.
// Fenwick stores:
//     Previous prefix-sum frequencies.
// Key observation:
//     pref[j]-pref[i]<k iff pref[i]>pref[j]-k.
// Extra technique:
//     Prefix inequality transformation.
// Returns:
//     Number of qualifying subarrays.
// Complexity:
//     O(n log n).
// Notes:
//     Implemented through the general inclusive range-sum counter safely at ll limits.
ll count_subarrays_sum_less_than_k(const vector<ll>&a, ll k) {
    vector<ll> prefix(1, 0);
    for(ll x:a) prefix.push_back(prefix.back() + x);
    vector<ll> values = prefix;
    sort(values.begin(), values.end());
    values.erase(unique(values.begin(), values.end()), values.end());
    FenwickTree fw(values.size());
    ll ans = 0, seen = 0;
    for(ll current:prefix) {
        i128 threshold = (i128) current - k;
        int lo = 0, hi = values.size();
        while(lo < hi) {
            int mid = (lo + hi) / 2;
            if((i128) values[mid] <= threshold) lo = mid + 1;
            else hi = mid;
        }
        ans += seen - fw.query(lo - 1);
        fw.update(lower_bound(values.begin(), values.end(), current) - values.begin(), 1);
        seen++;
    }
    return ans;
}

// Problem:
//     Count subarrays with positive sum.
// Fenwick stores:
//     Previous prefix-sum frequencies.
// Key observation:
//     pref[j]-pref[i]>0 iff pref[i]<pref[j].
// Extra technique:
//     Inversion-like counting of prefix sums.
// Returns:
//     Number of positive-sum subarrays.
// Complexity:
//     O(n log n).
// Notes:
//     Zero-sum subarrays are excluded by querying strictly smaller ranks.
ll count_positive_sum_subarrays(const vector<ll>&a) {
    vector<ll> prefix(1, 0);
    for(ll x:a) prefix.push_back(prefix.back() + x);
    auto [rank,values] = compress_values(prefix);
    FenwickTree fw(values.size());
    ll ans = 0;
    for(int p:rank) ans += fw.query(p - 1), fw.update(p, 1);
    return ans;
}

// Problem:
//     Count inversions in every length-k sliding window.
// Fenwick stores:
//     Frequencies of compressed values in the current window.
// Key observation:
//     Entering at right adds previous-greater; leaving at left removes later-smaller.
// Extra technique:
//     Sliding window with reversible additive frequencies.
// Returns:
//     n-k+1 inversion counts.
// Complexity:
//     O(n log n).
// Notes:
//     k must be between 1 and n.
vector<ll> sliding_window_inversions(const vector<ll>&a, int k) {
    assert(1 <= k && k <= (int)a.size());
    auto [rank,values] = compress_values(a);
    FenwickTree fw(values.size());
    vector<ll> ans;
    ll inversions = 0;
    for(int i = 0;i < (int) a.size();i++) {
        int inside = min(i, k - 1);
        inversions += inside - fw.query(rank[i]);
        fw.update(rank[i], 1);
        if(i >= k - 1) {
            ans.push_back(inversions);
            int leaving = i - k + 1;
            fw.update(rank[leaving], -1);
            inversions -= fw.query(rank[leaving] - 1);
        }
    }
    return ans;
}

// Problem:
//     For every length-k window compute sum |x-median|.
// Fenwick stores:
//     Counts and value sums by compressed coordinate.
// Key observation:
//     Cost splits into median*leftCount-leftSum and rightSum-median*rightCount.
// Extra technique:
//     Dynamic kth plus two Fenwick aggregates.
// Returns:
//     One cost per window, using the lower median.
// Complexity:
//     O(n log n).
// Notes:
//     Compression is possible because all array values are known beforehand.
vector<ll> sliding_window_median_cost(const vector<ll>&a, int k) {
    assert(1 <= k && k <= (int)a.size());
    FenwickOrderStatisticMultiset ms(a);
    vector<ll> ans;
    for(int i = 0;i < (int) a.size();i++) {
        ms.insert(a[i]);
        if(i >= k) ms.erase(a[i - k]);
        if(i >= k - 1) {
            ll med = ms.kth((k + 1) / 2), lc = ms.count_leq(med), ls = ms.sum_leq(med);
            ans.push_back(med * lc - ls + (ms.totalSum - ls) - med * (k - lc));
        }
    }
    return ans;
}

// Problem:
//     Count for each interval how many it contains and how many contain it.
// Fenwick stores:
//     Frequencies of right endpoints from already processed left endpoints.
// Key observation:
//     Intervals are 2D points (l,r); containment is dominance.
// Extra technique:
//     Two sweeps with duplicate groups.
// Returns:
//     {contains,containedBy}, excluding self but including identical other intervals.
// Complexity:
//     O(n log n).
// Notes:
//     Grouping identical intervals is required for symmetric duplicate answers.
pair<vector<ll>, vector<ll> > nested_interval_counts(const vector<pair<ll, ll> >&intervals) {
    int n = intervals.size();
    vector<ll> rights;
    for(auto [l,r]:intervals) rights.push_back(r);
    sort(rights.begin(), rights.end());
    rights.erase(unique(rights.begin(), rights.end()), rights.end());
    vector<int> order(n);
    iota(order.begin(), order.end(), 0);
    vector<ll> contains(n), contained(n);
    auto same = [&](int a, int b) { return intervals[a] == intervals[b]; };
    sort(order.begin(), order.end(), [&](int a, int b) {
        if(intervals[a].first != intervals[b].first) return intervals[a].first > intervals[b].first;
        return intervals[a].second < intervals[b].second;
    });
    FenwickTree fw(rights.size());
    for(int l = 0, r;l < n;l = r) {
        for(r = l;r < n && same(order[l], order[r]);r++);
        int copies = r - l;
        int p = lower_bound(rights.begin(), rights.end(), intervals[order[l]].second) - rights.begin();
        ll base = fw.query(p) + copies - 1;
        for(int i = l;i < r;i++) contains[order[i]] = base;
        fw.update(p, copies);
    }
    sort(order.begin(), order.end(), [&](int a, int b) {
        if(intervals[a].first != intervals[b].first) return intervals[a].first < intervals[b].first;
        return intervals[a].second > intervals[b].second;
    });
    fw = FenwickTree(rights.size());
    ll processed = 0;
    for(int l = 0, r;l < n;l = r) {
        for(r = l;r < n && same(order[l], order[r]);r++);
        int copies = r - l;
        int p = lower_bound(rights.begin(), rights.end(), intervals[order[l]].second) - rights.begin();
        ll base = processed - fw.query(p - 1) + copies - 1;
        for(int i = l;i < r;i++) contained[order[i]] = base;
        fw.update(p, copies);
        processed += copies;
    }
    return {contains, contained};
}

struct EulerTour {
    vector<int> tin, tout, parent, order;

    // Problem:
    //     Flatten every rooted subtree into one contiguous array interval.
    // Fenwick stores:
    //     Later structures index nodes by tin.
    // Key observation:
    //     DFS enters all descendants of u before leaving u.
    // Extra technique:
    //     Iterative Euler tour.
    // Returns:
    //     tin[u], tout[u] inclusive, parent and preorder.
    // Complexity:
    //     O(n).
    // Notes:
    //     Input is an undirected connected tree.
    EulerTour(const vector<vector<int> >&tree = {}, int root = 0) {
        int n = tree.size();
        tin.assign(n, -1);
        tout.assign(n, -1);
        parent.assign(n, -1);
        if(!n) return;
        vector<pair<int, int> > stack{{root, 0}};
        parent[root] = root;
        while(!stack.empty()) {
            auto&[u,next] = stack.back();
            if(!next && tin[u] == -1) tin[u] = order.size(), order.push_back(u);
            if(next == (int) tree[u].size()) {
                tout[u] = (int) order.size() - 1;
                stack.pop_back();
                continue;
            }
            int v = tree[u][next++];
            if(v == parent[u]) continue;
            if(parent[v] != -1) throw invalid_argument("EulerTour: graph is not a tree");
            parent[v] = u;
            stack.push_back({v, 0});
        }
        if((int) order.size() != n) throw invalid_argument("EulerTour: disconnected tree");
    }
};

struct SubtreeFenwick {
    EulerTour tour;
    vector<ll> value;
    FenwickTree fw;

    SubtreeFenwick(const vector<vector<int> >&tree, const vector<ll>&initial, int root = 0)
        : tour(tree, root), value(initial), fw((int) tree.size()) {
        assert(tree.size() == initial.size());
        for(int u = 0;u < (int) tree.size();u++) fw.update(tour.tin[u], value[u]);
    }

    // Problem:
    //     Assign a node value and query dynamic subtree sums.
    // Fenwick stores:
    //     Node values at Euler entry positions.
    // Key observation:
    //     Subtree(u) is [tin[u],tout[u]].
    // Extra technique:
    //     Euler tour flattening.
    // Returns:
    //     subtree_query returns the current subtree sum.
    // Complexity:
    //     O(log n) per update/query after O(n log n) initialization.
    // Notes:
    //     update_node is assignment, implemented through a delta.
    void update_node(int u, ll newValue) {
        fw.update(tour.tin[u], newValue - value[u]);
        value[u] = newValue;
    }

    ll subtree_query(int u) const { return fw.query(tour.tin[u], tour.tout[u]); }
};

struct SubtreeAddPointQuery {
    EulerTour tour;
    vector<ll> initial;
    RangeAddPointQueryFenwick fw;

    SubtreeAddPointQuery(const vector<vector<int> >&tree, const vector<ll>&initialValues, int root = 0)
        : tour(tree, root), initial(initialValues), fw(tree.size()) {
        assert(tree.size() == initialValues.size());
    }

    void update_subtree(int u, ll v) { fw.update(tour.tin[u], tour.tout[u], v); }
    ll query_node(int u) const { return initial[u] + fw.query(tour.tin[u]); }
};

struct SubtreeAddSubtreeSum {
    EulerTour tour;
    RangeFenwickTree fw;

    SubtreeAddSubtreeSum(const vector<vector<int> >&tree, const vector<ll>&initial, int root = 0)
        : tour(tree, root), fw(tree.size()) {
        assert(tree.size() == initial.size());
        for(int u = 0;u < (int) tree.size();u++) fw.update(tour.tin[u], tour.tin[u], initial[u]);
    }

    void update_subtree(int u, ll v) { fw.update(tour.tin[u], tour.tout[u], v); }
    ll query_subtree(int u) const { return fw.query(tour.tin[u], tour.tout[u]); }
};

struct HldFenwick {
    int n = 0, root = 0, timer = 0;
    vector<vector<int> > tree;
    vector<int> parent, depth, heavy, head, pos, size;
    vector<ll> value;
    FenwickTree fw;

    // Problem:
    //     Support point assignment and sum on any tree path.
    // Fenwick stores:
    //     Node values in heavy-light base-array order.
    // Key observation:
    //     Every path decomposes into O(log n) contiguous heavy-chain intervals.
    // Extra technique:
    //     Heavy-Light Decomposition.
    // Returns:
    //     path_query(u,v) returns the inclusive node-path sum.
    // Complexity:
    //     O(log n) update and O(log^2 n) path query.
    // Notes:
    //     Fenwick is suitable because the path operation is addition with an inverse.
    HldFenwick(const vector<vector<int> >&adjacency, const vector<ll>&initial, int rootNode = 0)
        : n(adjacency.size()), root(rootNode), tree(adjacency), parent(n, -1), depth(n), heavy(n, -1),
          head(n), pos(n), size(n, 1), value(initial), fw(n) {
        assert(adjacency.size() == initial.size());
        if(!n) return;
        vector<int> order{root};
        parent[root] = root;
        for(int at = 0;at < n;at++) {
            int u = order[at];
            for(int v:tree[u])
                if(v != parent[u]) {
                    if(parent[v] != -1) throw invalid_argument("HldFenwick: graph is not a tree");
                    parent[v] = u;
                    depth[v] = depth[u] + 1;
                    order.push_back(v);
                }
        }
        if((int) order.size() != n) throw invalid_argument("HldFenwick: disconnected tree");
        for(int at = n - 1;at > 0;at--) {
            int u = order[at], p = parent[u];
            size[p] += size[u];
            if(heavy[p] == -1 || size[u] > size[heavy[p]]) heavy[p] = u;
        }
        vector<pair<int, int> > stack{{root, root}};
        while(!stack.empty()) {
            auto [start,h] = stack.back();
            stack.pop_back();
            for(int u = start;u != -1;u = heavy[u]) {
                head[u] = h;
                pos[u] = timer++;
                for(int v:tree[u]) if(parent[v] == u && v != heavy[u]) stack.push_back({v, v});
            }
        }
        for(int u = 0;u < n;u++) fw.update(pos[u], value[u]);
    }

    void update_node(int u, ll newValue) {
        fw.update(pos[u], newValue - value[u]);
        value[u] = newValue;
    }

    ll path_query(int u, int v) const {
        ll ans = 0;
        while(head[u] != head[v]) {
            if(depth[head[u]] < depth[head[v]]) swap(u, v);
            ans += fw.query(pos[head[u]], pos[u]);
            u = parent[head[u]];
        }
        if(depth[u] > depth[v]) swap(u, v);
        return ans + fw.query(pos[u], pos[v]);
    }
};

// Problem:
//     For every node, count strict-greater ancestors on its root path.
// Fenwick stores:
//     Frequencies of compressed values currently on the DFS stack.
// Key observation:
//     Insert on entry and erase on exit, so the BIT represents exactly the ancestors.
// Extra technique:
//     DFS event stack and coordinate compression by value.
// Returns:
//     One ancestor-inversion count per node.
// Complexity:
//     O(n log n).
// Notes:
//     Here the Fenwick index represents VALUES, not Euler positions.
vector<ll> greater_ancestors(const vector<vector<int> >&tree,
                             const vector<ll>&value, int root = 0) {
    int n = tree.size();
    assert(value.size() == tree.size());
    if(!n) return {};
    EulerTour tour(tree, root);
    auto [rank,values] = compress_values(value);
    FenwickTree fw(values.size());
    vector<ll> ans(n);
    struct Event {
        int u, depth;
        bool leave;
    };
    vector<Event> stack{{root, 0, false}};
    while(!stack.empty()) {
        auto event = stack.back();
        stack.pop_back();
        int u = event.u;
        if(event.leave) {
            fw.update(rank[u], -1);
            continue;
        }
        ans[u] = event.depth - fw.query(rank[u]);
        fw.update(rank[u], 1);
        stack.push_back({u, event.depth, true});
        for(int i = (int) tree[u].size() - 1;i >= 0;i--) {
            int v = tree[u][i];
            if(tour.parent[v] == u) stack.push_back({v, event.depth + 1, false});
        }
    }
    return ans;
}

// Problem:
//     Convert a permutation to Lehmer code c[i]=number of later smaller values.
// Fenwick stores:
//     Frequencies of values already seen in a right-to-left sweep.
// Key observation:
//     c[i] is a strict prefix frequency query by value.
// Extra technique:
//     Permutation ranking.
// Returns:
//     Lehmer code with 0<=c[i]<n-i.
// Complexity:
//     O(n log n).
// Notes:
//     p must be a permutation of 0..n-1.
vector<int> lehmer_code(const vector<int>&p) {
    int n = p.size();
    FenwickTree fw(n);
    vector<char> seen(n);
    vector<int> code(n);
    for(int i = n - 1;i >= 0;i--) {
        assert(0 <= p[i] && p[i] < n && !seen[p[i]]);
        seen[p[i]] = 1;
        code[i] = fw.query(p[i] - 1);
        fw.update(p[i], 1);
    }
    return code;
}

// Problem:
//     Reconstruct a permutation from c[i]=number of later smaller values.
// Fenwick stores:
//     Availability frequency of each value 0..n-1.
// Key observation:
//     p[i] is the (c[i]+1)-th smallest currently unused value.
// Extra technique:
//     Fenwick kth/order statistics.
// Returns:
//     The unique permutation for a valid code.
// Complexity:
//     O(n log n).
// Notes:
//     This is the left-to-right Lehmer convention, not every inversion-sequence convention.
vector<int> permutation_from_lehmer(const vector<int>&code) {
    int n = code.size();
    FenwickTree fw(vector<ll>(n, 1));
    vector<int> p(n);
    for(int i = 0;i < n;i++) {
        assert(0 <= code[i] && code[i] < n-i);
        p[i] = fw.kth(code[i] + 1);
        fw.update(p[i], -1);
    }
    return p;
}

// Problem:
//     Kendall tau distance between two permutations of the same distinct values.
// Fenwick stores:
//     Used indirectly by inversion_count on mapped target positions.
// Key observation:
//     A disagreement is an inversion after mapping a-values to their positions in b.
// Extra technique:
//     Position mapping.
// Returns:
//     Number of pairwise order disagreements.
// Complexity:
//     O(n log n) expected.
// Notes:
//     Values must be unique and the sets equal.
ll kendall_tau(const vector<ll>&a, const vector<ll>&b) {
    assert(a.size() == b.size());
    unordered_map<ll, int> position;
    for(int i = 0;i < (int) b.size();i++)
        assert(!position.count(b[i])), position[b[i]] = i;
    vector<ll> mapped;
    for(ll x:a)
        assert(position.count(x)), mapped.push_back(position[x]);
    return inversion_count(mapped);
}

// Problem:
//     Minimum adjacent swaps needed to transform sequence a into sequence b.
// Fenwick stores:
//     Used by inversion_count on stable target occurrence positions.
// Key observation:
//     Match equal duplicates in occurrence order; crossings are required swaps.
// Extra technique:
//     Queues of target positions.
// Returns:
//     Minimum swaps, or -1 if the multisets differ.
// Complexity:
//     O(n log n) expected.
// Notes:
//     Stable duplicate matching is optimal.
ll minimum_adjacent_swaps(const vector<ll>&a, const vector<ll>&b) {
    if(a.size() != b.size()) return -1;
    unordered_map<ll, queue<int> > positions;
    for(int i = 0;i < (int) b.size();i++) positions[b[i]].push(i);
    vector<ll> mapped;
    for(ll x:a) {
        if(positions[x].empty()) return -1;
        mapped.push_back(positions[x].front());
        positions[x].pop();
    }
    return inversion_count(mapped);
}

ll minimum_adjacent_swaps(const string&a, const string&b) {
    vector<ll> x(a.begin(), a.end()), y(b.begin(), b.end());
    return minimum_adjacent_swaps(x, y);
}

// Problem:
//     Inversion counts after each left rotation of a permutation 0..n-1.
// Fenwick stores:
//     Only the initial inversion count.
// Key observation:
//     Moving first x to end removes x smaller-after pairs and adds n-1-x larger-before pairs.
// Extra technique:
//     O(1) mathematical update after Fenwick preprocessing.
// Returns:
//     Count before each of n rotations.
// Complexity:
//     O(n log n).
// Notes:
//     Formula is specific to a permutation of consecutive values 0..n-1.
vector<ll> rotation_inversion_counts(const vector<int>&p) {
    int n = p.size();
    vector<char> seen(n);
    vector<ll> asLl;
    for(int x:p)
        assert(0 <= x && x < n && !seen[x]), seen[x] = 1, asLl.push_back(x);
    vector<ll> ans(n);
    if(!n) return ans;
    ans[0] = inversion_count(asLl);
    for(int i = 1;i < n;i++) ans[i] = ans[i - 1] + (n - 1 - 2LL * p[i - 1]);
    return ans;
}

struct Point2D {
    ll x = 0, y = 0, weight = 1;
};

// Problem:
//     For each point count points with x'<=x and y'<=y, including itself.
// Fenwick stores:
//     Frequencies by compressed y among points with eligible x.
// Key observation:
//     Sort/sweep x; Fenwick handles the remaining y prefix.
// Extra technique:
//     Offline 2D dominance and equal-x batching.
// Returns:
//     Dominance counts in original order.
// Complexity:
//     O(n log n).
// Notes:
//     Insert an entire equal-x group before querying it because x equality is allowed.
vector<ll> dominance_2d(const vector<pair<ll, ll> >&points) {
    int n = points.size();
    vector<ll> ys;
    vector<int> order(n);
    iota(order.begin(), order.end(), 0);
    for(auto [x,y]:points) ys.push_back(y);
    sort(ys.begin(), ys.end());
    ys.erase(unique(ys.begin(), ys.end()), ys.end());
    sort(order.begin(), order.end(), [&](int a, int b) { return points[a] < points[b]; });
    FenwickTree fw(ys.size());
    vector<ll> ans(n);
    for(int l = 0, r;l < n;l = r) {
        for(r = l;r < n && points[order[r]].first == points[order[l]].first;r++);
        for(int i = l;i < r;i++) {
            int p = lower_bound(ys.begin(), ys.end(), points[order[i]].second) - ys.begin();
            fw.update(p, 1);
        }
        for(int i = l;i < r;i++) {
            int p = upper_bound(ys.begin(), ys.end(), points[order[i]].second) - ys.begin() - 1;
            ans[order[i]] = fw.query(p);
        }
    }
    return ans;
}

struct RectangleQuery {
    ll x1 = 0, y1 = 0, x2 = -1, y2 = -1;
};

// Problem:
//     Sum weights of static points inside each closed axis-aligned rectangle.
// Fenwick stores:
//     Weights by compressed y for points activated by the x sweep.
// Key observation:
//     Rectangle = prefix(x2,y-range)-prefix(x<x1,y-range).
// Extra technique:
//     Sweep-line events and inclusion-exclusion.
// Returns:
//     Weighted rectangle answers; weights=1 gives point counts.
// Complexity:
//     O((n+q) log(n+q)).
// Notes:
//     Strict-before and inclusive-after event types avoid x1-1 overflow.
vector<ll> offline_rectangle_sums(vector<Point2D> points,
                                  const vector<RectangleQuery>&queries) {
    struct Event {
        ll x, y1, y2;
        int id, sign;
        bool inclusive;
    };
    vector<ll> ys;
    for(auto p:points) ys.push_back(p.y);
    vector<Event> events;
    for(int i = 0;i < (int) queries.size();i++) {
        auto q = queries[i];
        if(q.x1 <= q.x2 && q.y1 <= q.y2) {
            events.push_back({q.x1, q.y1, q.y2, i, -1, false});
            events.push_back({q.x2, q.y1, q.y2, i, 1, true});
        }
    }
    sort(points.begin(), points.end(), [](Point2D a, Point2D b) { return a.x < b.x; });
    sort(ys.begin(), ys.end());
    ys.erase(unique(ys.begin(), ys.end()), ys.end());
    sort(events.begin(), events.end(), [](const Event&a, const Event&b) {
        if(a.x != b.x) return a.x < b.x;
        return a.inclusive < b.inclusive;
    });
    FenwickTree fw(ys.size());
    vector<ll> ans(queries.size());
    int at = 0;
    for(auto event:events) {
        while(at < (int) points.size() && (points[at].x < event.x ||
                                           (event.inclusive && points[at].x == event.x))) {
            int p = lower_bound(ys.begin(), ys.end(), points[at].y) - ys.begin();
            fw.update(p, points[at].weight);
            at++;
        }
        int l = lower_bound(ys.begin(), ys.end(), event.y1) - ys.begin();
        int r = upper_bound(ys.begin(), ys.end(), event.y2) - ys.begin() - 1;
        ans[event.id] += event.sign * fw.query(l, r);
    }
    return ans;
}

struct HorizontalSegment {
    ll x1 = 0, x2 = 0, y = 0;
};

struct VerticalSegment {
    ll x = 0, y1 = 0, y2 = 0;
};

// Problem:
//     Count intersections between closed horizontal and vertical segments.
// Fenwick stores:
//     Active horizontal segments by compressed y.
// Key observation:
//     Sweep x with add at x1, query at vertical x, remove after x2.
// Extra technique:
//     Sweep-line event ordering add < query < remove.
// Returns:
//     Number of boundary-touching intersections included.
// Complexity:
//     O((h+v) log(h+v)).
// Notes:
//     Changing same-x event order changes whether touching endpoints count.
ll orthogonal_intersections(const vector<HorizontalSegment>&horizontal,
                            const vector<VerticalSegment>&vertical) {
    struct Event {
        ll x, y1, y2;
        int type;
    };
    vector<Event> events;
    vector<ll> ys;
    for(auto s:horizontal) {
        if(s.x1 > s.x2) swap(s.x1, s.x2);
        ys.push_back(s.y);
        events.push_back({s.x1, s.y, s.y, 0});
        events.push_back({s.x2, s.y, s.y, 2});
    }
    for(auto s:vertical) {
        if(s.y1 > s.y2) swap(s.y1, s.y2);
        events.push_back({s.x, s.y1, s.y2, 1});
    }
    sort(ys.begin(), ys.end());
    ys.erase(unique(ys.begin(), ys.end()), ys.end());
    sort(events.begin(), events.end(), [](const Event&a, const Event&b) {
        if(a.x != b.x) return a.x < b.x;
        return a.type < b.type;
    });
    FenwickTree fw(ys.size());
    ll ans = 0;
    for(auto e:events) {
        if(e.type != 1) {
            int p = lower_bound(ys.begin(), ys.end(), e.y1) - ys.begin();
            fw.update(p, e.type == 0 ? 1 : -1);
        } else {
            int l = lower_bound(ys.begin(), ys.end(), e.y1) - ys.begin();
            int r = upper_bound(ys.begin(), ys.end(), e.y2) - ys.begin() - 1;
            ans += fw.query(l, r);
        }
    }
    return ans;
}

struct Point3D {
    ll x = 0, y = 0, z = 0;
};

// Problem:
//     For every point count points with x'<=x,y'<=y,z'<=z, including itself.
// Fenwick stores:
//     Total multiplicity by compressed z during CDQ merges sorted by y.
// Key observation:
//     Sorted order handles x; CDQ separates earlier/later points; BIT handles z.
// Extra technique:
//     Duplicate grouping, CDQ divide-and-conquer, additive rollback.
// Returns:
//     Dominance count in original order.
// Complexity:
//     O(n log^2 n) time and O(n) memory.
// Notes:
//     Every update made in a merge is removed afterward.
vector<ll> dominance_3d_cdq(const vector<Point3D>&points) {
    struct Node {
        ll x, y, z, count = 0, answer = 0;
        vector<int> ids;
    };
    int original = points.size();
    vector<int> order(original);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b) {
        return tie(points[a].x, points[a].y, points[a].z) < tie(points[b].x, points[b].y, points[b].z);
    });
    vector<Node> nodes;
    vector<ll> zs;
    for(int id:order) {
        if(nodes.empty() || tie(nodes.back().x, nodes.back().y, nodes.back().z)
           != tie(points[id].x, points[id].y, points[id].z))
            nodes.push_back({points[id].x, points[id].y, points[id].z, 0, 0, {}});
        nodes.back().count++;
        nodes.back().ids.push_back(id);
        zs.push_back(points[id].z);
    }
    sort(zs.begin(), zs.end());
    zs.erase(unique(zs.begin(), zs.end()), zs.end());
    FenwickTree fw(zs.size());
    vector<Node> buffer(nodes.size());
    function<void(int, int)> cdq = [&](int l, int r) {
        if(r - l <= 1) return;
        int mid = (l + r) / 2;
        cdq(l, mid);
        cdq(mid, r);
        int i = l;
        for(int j = mid;j < r;j++) {
            while(i < mid && nodes[i].y <= nodes[j].y) {
                int p = lower_bound(zs.begin(), zs.end(), nodes[i].z) - zs.begin();
                fw.update(p, nodes[i].count);
                i++;
            }
            int p = upper_bound(zs.begin(), zs.end(), nodes[j].z) - zs.begin() - 1;
            nodes[j].answer += fw.query(p);
        }
        for(int k = l;k < i;k++) {
            int p = lower_bound(zs.begin(), zs.end(), nodes[k].z) - zs.begin();
            fw.update(p, -nodes[k].count);
        }
        merge(nodes.begin() + l, nodes.begin() + mid, nodes.begin() + mid, nodes.begin() + r,
              buffer.begin() + l, [](const Node&a, const Node&b) { return a.y < b.y; });
        copy(buffer.begin() + l, buffer.begin() + r, nodes.begin() + l);
    };
    cdq(0, nodes.size());
    vector<ll> ans(original);
    for(auto&node:nodes) for(int id:node.ids) ans[id] = node.answer + node.count;
    return ans;
}

struct TimedRangeUpdate {
    int l = 0, r = -1;
    ll add = 0;
};

// Problem:
//     Each entity owns positions; find earliest update count when its accumulated sum reaches target.
// Fenwick stores:
//     Difference events for range updates applied up to the current tested time.
// Key observation:
//     With nonnegative additions, satisfaction is monotone in time.
// Extra technique:
//     Parallel binary search and range-add/point-query Fenwick.
// Returns:
//     Earliest number of applied events in [0,m], or -1 if never reached.
// Complexity:
//     O((m+total owned positions)*log n*log(m+1)).
// Notes:
//     Negative updates invalidate monotonicity. Positions may belong to multiple entities.
vector<int> earliest_owned_threshold(int positionCount,
                                     const vector<vector<int> >&ownedPositions, const vector<ll>&target,
                                     const vector<TimedRangeUpdate>&updates) {
    int q = ownedPositions.size(), m = updates.size();
    assert(target.size() == ownedPositions.size());
    vector<int> low(q, 0), high(q, m + 1);
    vector<char> done(q);
    for(int id = 0;id < q;id++) if(target[id] <= 0) high[id] = 0, done[id] = 1;
    while(true) {
        vector<vector<int> > bucket(m + 1);
        bool any = false;
        for(int id = 0;id < q;id++)
            if(!done[id] && high[id] - low[id] > 1) {
                bucket[(low[id] + high[id]) / 2].push_back(id);
                any = true;
            }
        if(!any) break;
        RangeAddPointQueryFenwick fw(positionCount);
        for(int time = 1;time <= m;time++) {
            auto event = updates[time - 1];
            assert(event.add >= 0);
            fw.update(event.l, event.r, event.add);
            for(int id:bucket[time]) {
                i128 sum = 0;
                for(int p:ownedPositions[id]) {
                    assert(0 <= p && p < positionCount);
                    sum += fw.query(p);
                }
                if(sum >= target[id]) high[id] = time;
                else low[id] = time;
            }
        }
    }
    vector<int> ans(q);
    for(int i = 0;i < q;i++) ans[i] = high[i] == m + 1 ? -1 : high[i];
    return ans;
}

/*
FUNCTION MAP
    FenwickTree                         point add + prefix/range sum + kth
    RangeAddPointQueryFenwick           range add + point value
    RangeFenwickTree                    range add + range sum
    PolynomialFenwickTree               constant/AP/quadratic range add + sum
    XorFenwickTree                      point xor + range xor
    FenwickMin / FenwickMax             monotone prefix min/max
    FenwickOrderStatisticMultiset       duplicates, kth, median, top-k sums
    FenwickTree2D                       dense point add + rectangle sum
    RangeFenwickTree2D                  dense rectangle add + rectangle sum
    SparseFenwickTree2D                 predeclared sparse points, rectangle sum

    inversion_count                     strict inversion pairs
    count_smaller_after                 suffix ranks per index
    count_greater_before                prefix ranks per index
    reverse_pairs                       a[i] > 2*a[j]
    count_increasing_triplets           a[i] < a[j] < a[k]
    count_crossings                     two-line edge crossings
    chord_intersections                 circular chord crossings
    josephus_order                      kth alive deletion
    remove_by_orders                    CSES-style list removals
    lis_*                               max/count DP indexed by value
    range_count_leq / range_sum_leq     offline threshold queries
    distinct_range_queries              last-occurrence right sweep
    exactly_once_range_queries          two-previous-occurrence events
    count_range_sum                     subarray-sum interval counting
    sliding_window_*                    inversion / median distance windows
    nested_interval_counts              contains / contained-by counts
    Subtree* / HldFenwick               Euler and heavy-light tree queries
    greater_ancestors                   DFS BIT indexed by values
    dominance_2d                        sort x + BIT y
    offline_rectangle_sums              four prefix sweep events
    orthogonal_intersections            add/query/remove sweep
    dominance_3d_cdq                    CDQ + BIT + rollback
    earliest_owned_threshold            parallel binary search + BIT

COMPLEXITY / MEMORY
    1D operation                        O(log n), O(n) memory
    linear build                        O(n)
    kth by BIT lifting                  O(log n)
    range variants                      O(log n), 2 or 4 arrays
    dense 2D operation                  O(log n log m), O(n*m) memory
    sparse 2D operation                 O(log^2 n), O(points log n) memory
    typical offline sweep               O((n+q) log(n+q))
    HLD path sum                        O(log^2 n)
    CDQ dominance                       O(n log^2 n)
    parallel binary search              O((updates+checks) log n log time)

STRICTNESS CHEAT SHEET
    count < x       lower_bound(x), query(position-1)
    count <= x      upper_bound(x), query(position-1)
    count > x       total-count<=x
    count >= x      total-count<x
    strict first sorted coordinate: query equal-key group, then insert group
    non-strict first coordinate:      insert equal-key group, then query group

WHEN ANOTHER TOOL IS BETTER
    static range sums                  prefix sums
    all offline range additions        difference array
    arbitrary range min/max            segment tree / sparse table
    static subarray kth                wavelet tree / persistent segment tree
    online unknown ordered keys        balanced BST or PBDS
    both query boundaries move freely  often Mo's algorithm
    complicated lazy range updates     segment tree

RECOGNITION EXERCISES (answer after the arrow)
 1. point add + range sum                         -> basic BIT, index=position, store=sum
 2. point assignment + sum                        -> basic BIT + mirror array delta
 3. range add + point value                       -> difference BIT
 4. range add + range sum                         -> two BITs
 5. add arithmetic progressions                   -> polynomial BIT
 6. kth alive index                               -> frequency BIT binary lifting
 7. static sums only                              -> prefix sums, no BIT
 8. arbitrary range minimum with assignments      -> segment tree, not BIT
 9. inversions with huge values                   -> compress values, store frequencies
10. [l,r] count values <=x                        -> sort x, BIT indexes positions
11. distinct values in [l,r]                      -> sweep r, store latest occurrences
12. subtree point updates and sums                -> Euler interval + basic BIT
13. subtree add and node value                     -> Euler + difference BIT
14. path sum with point changes                    -> HLD + BIT
15. ancestors with value <=x                       -> DFS stack, BIT indexes values
16. increasing-subsequence DP                     -> BIT indexes final value, stores max/sum
17. dynamic median on known values                 -> compressed frequency BIT
18. sum of largest k active scores                 -> count BIT + value-sum BIT
19. 2D dominance                                  -> sweep x, BIT indexes y
20. static rectangle point counts                 -> four x-prefix events + BIT y
21. horizontal/vertical intersections             -> add/query/remove x events + BIT y
22. 3D dominance                                  -> CDQ + BIT
23. earliest update reaching target                -> parallel binary search + BIT
24. static subarray kth                            -> wavelet/persistent tree, not plain BIT
25. online arbitrary keys without compression      -> balanced tree/PBDS

REAL PROBLEM CATALOG
    CSES Range Sum Queries II       -> FenwickTree + assign_value
    CSES List Removals              -> remove_by_orders / kth
    CSES Josephus Problem II        -> josephus_order
    CSES Distinct Values Queries    -> distinct_range_queries
    CSES Nested Ranges Count        -> nested_interval_counts
    CSES Sliding Window Cost        -> sliding_window_median_cost
    CSES Forest Queries II          -> FenwickTree2D
    SPOJ INVCNT                     -> inversion_count
    SPOJ DQUERY                     -> distinct_range_queries
    AtCoder DP Contest Q - Flowers  -> max_weight_increasing_subsequence

COMMON MISTAKES
    update means +=, not assignment; do not mix public 0-based with internal 1-based.
    query(-1)==0 and an r+1 update at n is intentionally safe in 1D wrappers.
    kth is invalid with negative frequencies or k outside [1,total].
    For arbitrary threshold queries choose lower_bound versus upper_bound carefully.
    Max/min BITs have no inverse: deletion and CDQ rollback do not work normally.
    Remove the old last occurrence in distinct queries and clean every CDQ update.
    Sweep event order decides whether boundary touching counts.
    Quadratic updates contain cubic expressions; switch coefficients to __int128 if needed.
*/

#ifndef FENWICK_TEMPLATE_NO_MAIN
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Example 1: point additions and range sums.
    FenwickTree fw(vector<ll>{3, 1, 4, 1, 5});
    fw.update(2, 6); // a[2] += 6
    cout << fw.query(1, 3) << '\n'; // 1+10+1 = 12

    // Example 2: range addition and range sum.
    RangeFenwickTree rangeFw(5);
    rangeFw.update(1, 3, 7);
    cout << rangeFw.query(2, 4) << '\n'; // 7+7+0 = 14

    // Example 3: inversion counting after compression.
    cout << inversion_count({3, 1, 2, 2}) << '\n';

    // Example 4: order statistics with duplicate values.
    FenwickOrderStatisticMultiset ms({10, 20, 30});
    ms.insert(20, 2);
    ms.insert(10);
    ms.insert(30);
    cout << ms.kth(3) << ' ' << ms.sum_k_largest(2) << '\n';

    // Example 5: an offline threshold query.
    vector<ThresholdRangeQuery> queries = {{1, 4, 3}, {0, 2, 10}};
    for(ll answer:range_count_leq({5, 1, 4, 2, 3}, queries)) cout << answer << ' ';
    cout << '\n';
}
#endif
