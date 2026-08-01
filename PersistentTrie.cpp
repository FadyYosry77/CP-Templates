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
        b = bits;
        clear();
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
        if(l > r)return 0;
 
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
        if(b < 63 && x == (1LL << b) - 1)return count(l, r);
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
};
