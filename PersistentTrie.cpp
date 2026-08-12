// PERSISTENT BINARY TRIE -- C++17
// B = configured number of bits, N = inserted values.
//
// # Build/rebuild from N values                 O(old nodes + N*B)
// # Nodes / space                               O(N*B), B+1 new nodes per insert
// # constructor / max_value / size              O(1)
// # reserve                                      O(current nodes) worst case
// # version_root / pop_version                  O(1), pop does not free nodes
// # clear                                        O(number of allocated nodes)
// # clone / count_nodes / count                  O(1)
// # insert_node / insert                         O(B)
// # count_equal / contains                       O(B)
// # count_less / count_less_equal                O(B)
// # count_greater / count_greater_equal          O(B)
// # count_between                                O(B)
// # kth_smallest / kth_largest                   O(B)
// # get_min / get_max                            O(B)
// # predecessor / floor_value                    O(B)
// # lower_bound_value / successor                O(B)
// # max_xor / min_xor                            O(B)
// # max_xor_element / min_xor_element            O(B)
// # max_xor_with_element / min_xor_with_element  O(B)
// # count_xor_less / count_xor_less_equal        O(B)
// # count_xor_greater / count_xor_greater_equal  O(B)
// # count_xor_between                            O(B)
// # kth_xor / kth_xor_element                    O(B)
// # kth_largest_xor / kth_largest_xor_element    O(B)
// # path_count / path_child_count                O(1)
// # path_count_equal                             O(B)
// # path_kth_smallest / path_kth_largest         O(B)
// # path_count_xor_less                          O(B)
// # path_max_xor / path_min_xor                  O(B)
// # maximum_subarray_xor                         O(N*B), O(N*B) extra space
// # count_subarrays_xor_less                     O(N*B), O(N*B) extra space
// # maximum_pair_xor                             O(N*B), O(N*B) extra space
// # count_pairs_xor_less                         O(N*B), O(N*B) extra space

/*
PERSISTENT BINARY TRIE RECOGNITION GUIDE

HOW TO RECOGNIZE A PERSISTENT-TRIE PROBLEM

1. The input is an array and every query chooses a subarray [l,r].
   The query asks about the VALUES in that range, especially:
   - maximum/minimum x^a[i];
   - how many (x^a[i]) are below k;
   - kth smallest value or kth smallest XOR;
   - frequency, predecessor, successor, or numeric rank.

   Build root[i] from the first i values.  The multiset of [l,r] is the
   difference root[r+1]-root[l], exactly like prefix sums but at every bit.

2. The statement asks to keep or query old versions after insertions.
   Every insertion changes only one root-to-leaf binary path, so path copying
   creates a new version in O(bits) without changing older versions.

3. Values live naturally in binary and n,q are large enough that scanning a
   range is impossible.  Typical constraints are n,q around 2e5 and values
   below 2^30/2^60, making O(bits) per query appropriate.

4. Values are placed on tree vertices and the query is on a path u..v.
   Build one version for every root-to-vertex path.  With w=lca(u,v):
       path multiset = root[u]+root[v]-root[w]-root[parent[w]].
   This gives path kth/frequency/XOR queries in O(bits), after LCA work.

5. Prefix XOR appears in the transformation.
   Since xor(l..r)=prefix[r]^prefix[l-1], inserting prefix XORs into a trie
   often turns maximum/counting subarray-XOR problems into ordinary queries.

WHEN NOT TO USE IT
- Only one changing set and no old-version/range query: use a normal bit trie.
- Queries are sums/minimums by array index, not binary value decisions: use a
  Fenwick/segment tree or persistent segment tree.
- Values change online at arbitrary indices: prefix versions no longer model
  the array directly; consider a trie of segment trees or an offline method.

CONTAINS
- One immutable prefix version per inserted array value.
- Range frequency, order statistics, predecessor/successor and bounds.
- Minimum/maximum/kth XOR values and the elements producing them.
- XOR threshold and interval counting.
- Four-root tree-path queries for kth, frequency, min XOR and max XOR.
- Standard maximum/counting recipes for subarray XOR and pair XOR.

CONVENTIONS
- Stored values are nonnegative and fit in b bits; 1<=b<=63.
- Array ranges [l,r] are zero-indexed and inclusive.
- k in kth queries is one-based.
- -1 means an invalid/empty answer because valid stored values are nonnegative.
- Each insert creates b+1 nodes: memory is O(numberOfValues*b).
- Prefix range queries require roots made by consecutive insert().
- Tree versions are made with insert_node(parentRoot,value).

COMMON CHOICES
- Values below 2^20: use b=20.
- int values up to 1e9: use b=31.
- Nonnegative signed long long values: use b=63.
*/

#include <bits/stdc++.h>
using namespace std;

using ll = long long;
#define pb push_back

struct persistentTrie {
    struct node {
        int to[2], cnt;
 
        node() {
            to[0] = to[1] = cnt = 0;
        }
    };
 
    int b;
    vector<node> tr;
    vector<int> root;
 
    persistentTrie(int bits = 31) {
        assert(1 <= bits && bits <= 63);
        b = bits;
        clear();
    }

    // Largest representable nonnegative value. O(1)
    ll max_value() const {
        return b == 63 ? LLONG_MAX : (1LL << b) - 1;
    }

    // Reserves enough memory for expected inserted values. O(1)
    void reserve(int values) {
        assert(values >= 0);
        tr.reserve(1ULL + (unsigned long long)values * (b + 1));
        root.reserve(values + 1);
    }
 
    // Removes all values and creates an empty version.
    void clear() {
        tr.clear();
        root.clear();
        tr.pb(node());
        root.pb(0);
    }
 
    // Creates a copy of an existing node.
    int clone(int u) {
        tr.pb(tr[u]);
        return tr.size() - 1;
    }
 
    // Returns how many values exist between two trie nodes.
    int count_nodes(int l,int r) const {
        return tr[r].cnt - tr[l].cnt;
    }
 
    // Creates a new root after inserting x into an old version.
    int insert_node(int old,ll x) {
        assert(0 <= x && x <= max_value());
        assert(0 <= old && old < (int)tr.size());
        int rt = clone(old);
        int u = rt, v = old;
        tr[u].cnt++;
 
        for(int i = b - 1;i >= 0;i--) {
            int bit = x >> i & 1;
            int nxt = clone(tr[v].to[bit]);
 
            tr[u].to[bit] = nxt;
            u = nxt;
            v = tr[v].to[bit];
            tr[u].cnt++;
        }
 
        return rt;
    }
 
    // Inserts x and creates the next version.
    void insert(ll x) {
        root.pb(insert_node(root.back(), x));
    }
 
    // Builds all persistent versions from an array.
    void build(const vector<ll>&a) {
        clear();
        for(ll x:a)insert(x);
    }
 
    // Returns the number of inserted values.
    int size() const {
        return root.size() - 1;
    }
 
    // Returns the number of values in index range [l,r].
    int count(int l,int r) const {
        if(l > r)return 0;
        return tr[root[r + 1]].cnt - tr[root[l]].cnt;
    }
 
    // Returns how many times x appears in index range [l,r].
    int count_equal(int l,int r,ll x) const {
        if(l > r || x < 0 || x > max_value())return 0;
 
        int u = root[l], v = root[r + 1];
 
        for(int i = b - 1;i >= 0;i--) {
            int bit = x >> i & 1;
            u = tr[u].to[bit];
            v = tr[v].to[bit];
        }
 
        return tr[v].cnt - tr[u].cnt;
    }
 
    // Checks whether x exists in index range [l,r].
    bool contains(int l,int r,ll x) const {
        return count_equal(l, r, x) > 0;
    }
 
    // Returns the maximum possible value of x^a[i] for l<=i<=r.
    ll max_xor(int l,int r,ll x) const {
        if(l > r)return -1;
 
        int u = root[l], v = root[r + 1];
        ll ans = 0;
 
        for(int i = b - 1;i >= 0;i--) {
            int bit = x >> i & 1;
            int want = bit ^ 1;
            int lu = tr[u].to[want];
            int rv = tr[v].to[want];
 
            if(tr[rv].cnt - tr[lu].cnt) {
                ans |= 1LL << i;
                u = lu;
                v = rv;
            } else {
                u = tr[u].to[bit];
                v = tr[v].to[bit];
            }
        }
 
        return ans;
    }
 
    // Returns the array value that produces the maximum XOR with x.
    ll max_xor_element(int l,int r,ll x) const {
        if(l > r)return -1;
 
        int u = root[l], v = root[r + 1];
        ll ans = 0;
 
        for(int i = b - 1;i >= 0;i--) {
            int bit = x >> i & 1;
            int want = bit ^ 1;
            int lu = tr[u].to[want];
            int rv = tr[v].to[want];
            int take;
 
            if(tr[rv].cnt - tr[lu].cnt)take = want;
            else take = bit;
 
            if(take)ans |= 1LL << i;
 
            u = tr[u].to[take];
            v = tr[v].to[take];
        }
 
        return ans;
    }
 
    // Returns the minimum possible value of x^a[i] for l<=i<=r.
    ll min_xor(int l,int r,ll x) const {
        if(l > r)return -1;
 
        int u = root[l], v = root[r + 1];
        ll ans = 0;
 
        for(int i = b - 1;i >= 0;i--) {
            int bit = x >> i & 1;
            int lu = tr[u].to[bit];
            int rv = tr[v].to[bit];
 
            if(tr[rv].cnt - tr[lu].cnt) {
                u = lu;
                v = rv;
            } else {
                ans |= 1LL << i;
                u = tr[u].to[bit ^ 1];
                v = tr[v].to[bit ^ 1];
            }
        }
 
        return ans;
    }
 
    // Returns the array value that produces the minimum XOR with x.
    ll min_xor_element(int l,int r,ll x) const {
        if(l > r)return -1;
 
        int u = root[l], v = root[r + 1];
        ll ans = 0;
 
        for(int i = b - 1;i >= 0;i--) {
            int bit = x >> i & 1;
            int lu = tr[u].to[bit];
            int rv = tr[v].to[bit];
            int take;
 
            if(tr[rv].cnt - tr[lu].cnt)take = bit;
            else take = bit ^ 1;
 
            if(take)ans |= 1LL << i;
 
            u = tr[u].to[take];
            v = tr[v].to[take];
        }
 
        return ans;
    }
 
    // Counts values where (a[i]^x)<k in index range [l,r].
    int count_xor_less(int l,int r,ll x,ll k) const {
        if(l > r || k <= 0)return 0;
        if(b < 63 && k >= (1LL << b))return count(l, r);
 
        int u = root[l], v = root[r + 1];
        int ans = 0;
 
        for(int i = b - 1;i >= 0;i--) {
            int xb = x >> i & 1;
            int kb = k >> i & 1;
 
            if(kb) {
                int lu = tr[u].to[xb];
                int rv = tr[v].to[xb];
 
                ans += tr[rv].cnt - tr[lu].cnt;
 
                u = tr[u].to[xb ^ 1];
                v = tr[v].to[xb ^ 1];
            } else {
                u = tr[u].to[xb];
                v = tr[v].to[xb];
            }
        }
 
        return ans;
    }
 
    // Counts values smaller than x in index range [l,r].
    int count_less(int l,int r,ll x) const {
        return count_xor_less(l, r, 0, x);
    }
 
    // Counts values smaller than or equal to x in [l,r].
    int count_less_equal(int l,int r,ll x) const {
        if(x < 0)return 0;
        if(x >= max_value())return count(l, r);
        return count_less(l, r, x + 1);
    }
 
    // Returns the kth smallest value in [l,r], where k is one-based.
    ll kth_smallest(int l,int r,int k) const {
        if(k < 1 || k > count(l, r))return -1;
 
        int u = root[l], v = root[r + 1];
        ll ans = 0;
 
        for(int i = b - 1;i >= 0;i--) {
            int lu = tr[u].to[0];
            int rv = tr[v].to[0];
            int left = tr[rv].cnt - tr[lu].cnt;
 
            if(k <= left) {
                u = lu;
                v = rv;
            } else {
                k -= left;
                ans |= 1LL << i;
                u = tr[u].to[1];
                v = tr[v].to[1];
            }
        }
 
        return ans;
    }
 
    // Returns the kth largest value in [l,r], where k is one-based.
    ll kth_largest(int l,int r,int k) const {
        int sz = count(l, r);
        if(k < 1 || k > sz)return -1;
        return kth_smallest(l, r, sz - k + 1);
    }
 
    // Returns the kth smallest XOR result among x^a[i] in [l,r].
    ll kth_xor(int l,int r,ll x,int k) const {
        if(k < 1 || k > count(l, r))return -1;
 
        int u = root[l], v = root[r + 1];
        ll ans = 0;
 
        for(int i = b - 1;i >= 0;i--) {
            int bit = x >> i & 1;
            int lu = tr[u].to[bit];
            int rv = tr[v].to[bit];
            int same = tr[rv].cnt - tr[lu].cnt;
 
            if(k <= same) {
                u = lu;
                v = rv;
            } else {
                k -= same;
                ans |= 1LL << i;
                u = tr[u].to[bit ^ 1];
                v = tr[v].to[bit ^ 1];
            }
        }
 
        return ans;
    }
 
    // Returns the minimum stored value in [l,r].
    ll get_min(int l,int r) const {
        return kth_smallest(l, r, 1);
    }
 
    // Returns the maximum stored value in [l,r].
    ll get_max(int l,int r) const {
        return kth_largest(l, r, 1);
    }


    // Removes the newest prefix version. Nodes remain allocated. O(1)
    bool pop_version() {
        if(root.size() == 1)return false;
        root.pop_back();
        return true;
    }

    // Returns the trie-node root of prefix version v. O(1)
    int version_root(int v) const {
        assert(0 <= v && v < (int)root.size());
        return root[v];
    }

    // Counts values in numeric interval [low,high]. O(b)
    int count_between(int l,int r,ll low,ll high) const {
        if(l > r || low > high)return 0;
        return count_less_equal(l, r, high) - count_less(l, r, low);
    }

    // Counts values strictly greater than x. O(b)
    int count_greater(int l,int r,ll x) const {
        return count(l, r) - count_less_equal(l, r, x);
    }

    // Counts values greater than or equal to x. O(b)
    int count_greater_equal(int l,int r,ll x) const {
        return count(l, r) - count_less(l, r, x);
    }

    // Largest value strictly smaller than x, or -1. O(b)
    ll predecessor(int l,int r,ll x) const {
        int rank = count_less(l, r, x);
        return rank ? kth_smallest(l, r, rank) : -1;
    }

    // Largest value <=x, or -1. O(b)
    ll floor_value(int l,int r,ll x) const {
        int rank = count_less_equal(l, r, x);
        return rank ? kth_smallest(l, r, rank) : -1;
    }

    // Smallest value >=x, or -1. O(b)
    ll lower_bound_value(int l,int r,ll x) const {
        int rank = count_less(l, r, x) + 1;
        return rank <= count(l, r) ? kth_smallest(l, r, rank) : -1;
    }

    // Smallest value strictly greater than x, or -1. O(b)
    ll successor(int l,int r,ll x) const {
        int rank = count_less_equal(l, r, x) + 1;
        return rank <= count(l, r) ? kth_smallest(l, r, rank) : -1;
    }

    // Counts values where (a[i]^x)<=k. O(b)
    int count_xor_less_equal(int l,int r,ll x,ll k) const {
        if(k < 0)return 0;
        if(k >= max_value())return count(l, r);
        return count_xor_less(l, r, x, k + 1);
    }

    // Counts values where (a[i]^x)>k. O(b)
    int count_xor_greater(int l,int r,ll x,ll k) const {
        return count(l, r) - count_xor_less_equal(l, r, x, k);
    }

    // Counts values where (a[i]^x)>=k. O(b)
    int count_xor_greater_equal(int l,int r,ll x,ll k) const {
        return count(l, r) - count_xor_less(l, r, x, k);
    }

    // Counts values where low<=(a[i]^x)<=high. O(b)
    int count_xor_between(int l,int r,ll x,ll low,ll high) const {
        if(l > r || low > high)return 0;
        return count_xor_less_equal(l, r, x, high)
             - count_xor_less(l, r, x, low);
    }

    // Returns the element producing the kth smallest XOR with x. O(b)
    ll kth_xor_element(int l,int r,ll x,int k) const {
        if(k < 1 || k > count(l, r))return -1;

        int u = root[l], v = root[r + 1];
        ll element = 0;

        for(int i = b - 1;i >= 0;i--) {
            int bit = x >> i & 1;
            int lu = tr[u].to[bit];
            int rv = tr[v].to[bit];
            int same = count_nodes(lu, rv);
            int take;

            if(k <= same)take = bit;
            else {
                k -= same;
                take = bit ^ 1;
            }

            if(take)element |= 1LL << i;
            u = tr[u].to[take];
            v = tr[v].to[take];
        }

        return element;
    }

    // Returns the kth largest XOR result. O(b)
    ll kth_largest_xor(int l,int r,ll x,int k) const {
        int sz = count(l, r);
        return k < 1 || k > sz ? -1 : kth_xor(l, r, x, sz - k + 1);
    }

    // Returns the element producing the kth largest XOR. O(b)
    ll kth_largest_xor_element(int l,int r,ll x,int k) const {
        int sz = count(l, r);
        return k < 1 || k > sz ? -1
                               : kth_xor_element(l, r, x, sz - k + 1);
    }

    // Returns {maximum XOR, element producing it}. O(b)
    pair<ll,ll> max_xor_with_element(int l,int r,ll x) const {
        if(l > r)return {-1,-1};
        return {max_xor(l,r,x),max_xor_element(l,r,x)};
    }

    // Returns {minimum XOR, element producing it}. O(b)
    pair<ll,ll> min_xor_with_element(int l,int r,ll x) const {
        if(l > r)return {-1,-1};
        return {min_xor(l,r,x),min_xor_element(l,r,x)};
    }

    /*
    TREE-PATH MODE

    Build one trie root per tree vertex:
        vertexRoot[v] = insert_node(vertexRoot[parent[v]], value[v]);

    If w=lca(u,v), the path multiset is:
        root[u] + root[v] - root[w] - root[parent[w]]

    Pass 0 for root[parent[w]] when w is the tree root.
    */

    // Number of values represented by four inclusion-exclusion roots. O(1)
    int path_count(int ru,int rv,int rw,int rpw) const {
        return tr[ru].cnt + tr[rv].cnt - tr[rw].cnt - tr[rpw].cnt;
    }

    // Count contributed by one child bit under four path roots. O(1)
    int path_child_count(int ru,int rv,int rw,int rpw,int bit) const {
        return path_count(tr[ru].to[bit], tr[rv].to[bit],
                          tr[rw].to[bit], tr[rpw].to[bit]);
    }

    // Counts occurrences of x on a tree path. O(b)
    int path_count_equal(int ru,int rv,int rw,int rpw,ll x) const {
        for(int i = b - 1;i >= 0;i--) {
            int bit = x >> i & 1;
            ru = tr[ru].to[bit];
            rv = tr[rv].to[bit];
            rw = tr[rw].to[bit];
            rpw = tr[rpw].to[bit];
        }
        return path_count(ru,rv,rw,rpw);
    }

    // kth smallest value on a tree path, one-based. O(b)
    ll path_kth_smallest(int ru,int rv,int rw,int rpw,int k) const {
        if(k < 1 || k > path_count(ru,rv,rw,rpw))return -1;
        ll ans = 0;

        for(int i = b - 1;i >= 0;i--) {
            int left = path_child_count(ru,rv,rw,rpw,0);
            int take = 0;

            if(k > left) {
                k -= left;
                take = 1;
                ans |= 1LL << i;
            }

            ru = tr[ru].to[take];
            rv = tr[rv].to[take];
            rw = tr[rw].to[take];
            rpw = tr[rpw].to[take];
        }

        return ans;
    }

    // kth largest value on a tree path, one-based. O(b)
    ll path_kth_largest(int ru,int rv,int rw,int rpw,int k) const {
        int sz = path_count(ru,rv,rw,rpw);
        return k < 1 || k > sz ? -1
                               : path_kth_smallest(ru,rv,rw,rpw,sz-k+1);
    }

    // Counts path values y where (y^x)<k. O(b)
    int path_count_xor_less(int ru,int rv,int rw,int rpw,ll x,ll k) const {
        if(k <= 0)return 0;
        if(k >= max_value()) {
            if(k > max_value())return path_count(ru,rv,rw,rpw);
        }

        int ans = 0;

        for(int i = b - 1;i >= 0;i--) {
            int xb = x >> i & 1;
            int kb = k >> i & 1;

            if(kb) {
                ans += path_child_count(ru,rv,rw,rpw,xb);
                ru = tr[ru].to[xb ^ 1];
                rv = tr[rv].to[xb ^ 1];
                rw = tr[rw].to[xb ^ 1];
                rpw = tr[rpw].to[xb ^ 1];
            } else {
                ru = tr[ru].to[xb];
                rv = tr[rv].to[xb];
                rw = tr[rw].to[xb];
                rpw = tr[rpw].to[xb];
            }
        }

        return ans;
    }

    // Maximum x^value on a tree path. O(b)
    ll path_max_xor(int ru,int rv,int rw,int rpw,ll x) const {
        if(!path_count(ru,rv,rw,rpw))return -1;
        ll ans = 0;

        for(int i = b - 1;i >= 0;i--) {
            int bit = x >> i & 1;
            int take = bit ^ 1;

            if(!path_child_count(ru,rv,rw,rpw,take))take = bit;
            else ans |= 1LL << i;

            ru = tr[ru].to[take];
            rv = tr[rv].to[take];
            rw = tr[rw].to[take];
            rpw = tr[rpw].to[take];
        }

        return ans;
    }

    // Minimum x^value on a tree path. O(b)
    ll path_min_xor(int ru,int rv,int rw,int rpw,ll x) const {
        if(!path_count(ru,rv,rw,rpw))return -1;
        ll ans = 0;

        for(int i = b - 1;i >= 0;i--) {
            int bit = x >> i & 1;
            int take = bit;

            if(!path_child_count(ru,rv,rw,rpw,take)) {
                take ^= 1;
                ans |= 1LL << i;
            }

            ru = tr[ru].to[take];
            rv = tr[rv].to[take];
            rw = tr[rw].to[take];
            rpw = tr[rpw].to[take];
        }

        return ans;
    }

    // Maximum XOR of a nonempty subarray. O(n*b)
    static ll maximum_subarray_xor(const vector<ll>&a,int bits = 63) {
        if(a.empty())return 0;
        persistentTrie trie(bits);
        trie.reserve(a.size() + 1);
        trie.insert(0);

        ll prefix = 0, ans = 0;

        for(ll x:a) {
            prefix ^= x;
            ans = max(ans,trie.max_xor(0,trie.size()-1,prefix));
            trie.insert(prefix);
        }

        return ans;
    }

    // Counts subarrays with XOR<k. O(n*b)
    static long long count_subarrays_xor_less(
            const vector<ll>&a,ll k,int bits = 63) {
        persistentTrie trie(bits);
        trie.reserve(a.size() + 1);
        trie.insert(0);

        ll prefix = 0;
        long long ans = 0;

        for(ll x:a) {
            prefix ^= x;
            ans += trie.count_xor_less(0,trie.size()-1,prefix,k);
            trie.insert(prefix);
        }

        return ans;
    }

    // Maximum XOR of two different array positions, or -1 for size<2. O(n*b)
    static ll maximum_pair_xor(const vector<ll>&a,int bits = 63) {
        if(a.size() < 2)return -1;
        persistentTrie trie(bits);
        trie.reserve(a.size());
        trie.insert(a[0]);

        ll ans = 0;

        for(int i = 1;i < (int)a.size();i++) {
            ans = max(ans,trie.max_xor(0,trie.size()-1,a[i]));
            trie.insert(a[i]);
        }

        return ans;
    }

    // Counts index pairs i<j with (a[i]^a[j])<k. O(n*b)
    static long long count_pairs_xor_less(
            const vector<ll>&a,ll k,int bits = 63) {
        persistentTrie trie(bits);
        trie.reserve(a.size());

        long long ans = 0;

        for(ll x:a) {
            if(trie.size())
                ans += trie.count_xor_less(0,trie.size()-1,x,k);
            trie.insert(x);
        }

        return ans;
    }

};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    vector<ll> a = {5,1,7,4,3};
    persistentTrie trie(4);
    trie.reserve(a.size());
    trie.build(a);

    // Array range [1,4] contains {1,7,4,3}.
    assert(trie.count(1,4) == 4);
    assert(trie.count_equal(1,4,7) == 1);
    assert(trie.kth_smallest(1,4,2) == 3);
    assert(trie.count_between(1,4,3,5) == 2);
    assert(trie.predecessor(1,4,4) == 3);
    assert(trie.lower_bound_value(1,4,4) == 4);

    auto [bestXor,bestElement] = trie.max_xor_with_element(1,4,2);
    assert(bestXor == 6 && bestElement == 4);
    assert(trie.min_xor(1,4,2) == 1);
    assert(trie.count_xor_less(1,4,2,5) == 2);
    assert(trie.kth_xor(1,4,2,3) == 5);
    assert(trie.kth_xor_element(1,4,2,3) == 7);

    cout << "range max xor: " << bestXor << '\n';
    cout << "2nd smallest: " << trie.kth_smallest(1,4,2) << '\n';

    // Tree: 0(5), children 1(1),2(7), and 3(4) is child of 1.
    vector<int> treeRoot(4);
    treeRoot[0] = trie.insert_node(0,5);
    treeRoot[1] = trie.insert_node(treeRoot[0],1);
    treeRoot[2] = trie.insert_node(treeRoot[0],7);
    treeRoot[3] = trie.insert_node(treeRoot[1],4);

    // Path 3..2 has values {4,1,5,7}; lca=0, parent(lca)=none/root 0.
    int ru = treeRoot[3], rv = treeRoot[2], rw = treeRoot[0], rpw = 0;
    assert(trie.path_count(ru,rv,rw,rpw) == 4);
    assert(trie.path_kth_smallest(ru,rv,rw,rpw,2) == 4);
    assert(trie.path_max_xor(ru,rv,rw,rpw,2) == 7);
    assert(trie.path_count_xor_less(ru,rv,rw,rpw,2,5) == 1);

    cout << "tree path max xor: "
         << trie.path_max_xor(ru,rv,rw,rpw,2) << '\n';

    cout << "maximum subarray xor: "
         << persistentTrie::maximum_subarray_xor(a,4) << '\n';
    cout << "subarrays with xor < 4: "
         << persistentTrie::count_subarrays_xor_less(a,4,4) << '\n';
}
