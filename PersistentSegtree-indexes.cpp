#include<bits/stdc++.h>
using namespace std;

#define ll long long
#define pb push_back
#define all(v) (v).begin(),(v).end()

const ll inf = 4e18;

//====================================================
// Persistent Segment Tree - Index Version
// 0-indexed range: [0,n-1]
//====================================================

class PersistentSegtree {
private:
    struct Node {
        ll sum = 0;
        int l = 0, r = 0;
    };

    int n, timer;
    vector<Node> tree;

    Node join(int l, int r) {
        return Node{tree[l].sum + tree[r].sum, l, r};
    }

    int build_empty(int tl, int tr) {
        int cur = timer++;

        if(tl == tr) {
            tree[cur] = {0, 0, 0};
            return cur;
        }

        int mid = (tl + tr) / 2;
        int left = build_empty(tl, mid);
        int right = build_empty(mid + 1, tr);

        tree[cur] = join(left, right);

        return cur;
    }

    int build(int tl, int tr, const vector<ll>&arr) {
        int cur = timer++;

        if(tl == tr) {
            tree[cur] = {arr[tl], 0, 0};
            return cur;
        }

        int mid = (tl + tr) / 2;
        int left = build(tl, mid, arr);
        int right = build(mid + 1, tr, arr);

        tree[cur] = join(left, right);

        return cur;
    }

    int add(int v, int pos,ll val, int tl, int tr) {
        int cur = timer++;

        if(tl == tr) {
            tree[cur] = {tree[v].sum + val, 0, 0};
            return cur;
        }

        int mid = (tl + tr) / 2;

        if(pos <= mid) {
            int left = add(tree[v].l, pos, val, tl, mid);
            tree[cur] = join(left, tree[v].r);
        } else {
            int right = add(tree[v].r, pos, val, mid + 1, tr);
            tree[cur] = join(tree[v].l, right);
        }

        return cur;
    }

    int set_value(int v, int pos,ll val, int tl, int tr) {
        int cur = timer++;

        if(tl == tr) {
            tree[cur] = {val, 0, 0};
            return cur;
        }

        int mid = (tl + tr) / 2;

        if(pos <= mid) {
            int left = set_value(tree[v].l, pos, val, tl, mid);
            tree[cur] = join(left, tree[v].r);
        } else {
            int right = set_value(tree[v].r, pos, val, mid + 1, tr);
            tree[cur] = join(tree[v].l, right);
        }

        return cur;
    }

    ll query(int v, int ql, int qr, int tl, int tr) {
        if(!v)return 0;
        if(qr < tl || tr < ql)return 0;
        if(ql <= tl && tr <= qr)return tree[v].sum;

        int mid = (tl + tr) / 2;

        return query(tree[v].l, ql, qr, tl, mid)
               + query(tree[v].r, ql, qr, mid + 1, tr);
    }

    int kth_single(int v,ll k, int tl, int tr) {
        if(tl == tr)return tl;

        int mid = (tl + tr) / 2;
        ll left = tree[tree[v].l].sum;

        if(k <= left)return kth_single(tree[v].l, k, tl, mid);
        return kth_single(tree[v].r, k - left, mid + 1, tr);
    }

    int kth_range(int rroot, int lroot,ll k, int tl, int tr) {
        if(tl == tr)return tl;

        int mid = (tl + tr) / 2;
        ll left = tree[tree[rroot].l].sum - tree[tree[lroot].l].sum;

        if(k <= left) {
            return kth_range(tree[rroot].l, tree[lroot].l, k, tl, mid);
        } else {
            return kth_range(tree[rroot].r, tree[lroot].r, k - left, mid + 1, tr);
        }
    }

    int kth_path(int uroot, int vroot, int lcroot, int plcroot,ll k, int tl, int tr) {
        if(tl == tr)return tl;

        int mid = (tl + tr) / 2;

        ll left =
                tree[tree[uroot].l].sum +
                tree[tree[vroot].l].sum -
                tree[tree[lcroot].l].sum -
                tree[tree[plcroot].l].sum;

        if(k <= left) {
            return kth_path(
                tree[uroot].l,
                tree[vroot].l,
                tree[lcroot].l,
                tree[plcroot].l,
                k,
                tl,
                mid
            );
        } else {
            return kth_path(
                tree[uroot].r,
                tree[vroot].r,
                tree[lcroot].r,
                tree[plcroot].r,
                k - left,
                mid + 1,
                tr
            );
        }
    }

    ll path_query(int uroot, int vroot, int lcroot, int plcroot, int ql, int qr, int tl, int tr) {
        if(qr < tl || tr < ql)return 0;

        if(ql <= tl && tr <= qr) {
            return tree[uroot].sum +
                   tree[vroot].sum -
                   tree[lcroot].sum -
                   tree[plcroot].sum;
        }

        int mid = (tl + tr) / 2;

        return path_query(
                   tree[uroot].l,
                   tree[vroot].l,
                   tree[lcroot].l,
                   tree[plcroot].l,
                   ql,
                   qr,
                   tl,
                   mid
               )
               + path_query(
                   tree[uroot].r,
                   tree[vroot].r,
                   tree[lcroot].r,
                   tree[plcroot].r,
                   ql,
                   qr,
                   mid + 1,
                   tr
               );
    }

    int first_active(int v, int ql, int qr, int tl, int tr) {
        if(!v || tree[v].sum <= 0 || qr < tl || tr < ql)return -1;
        if(tl == tr)return tl;

        int mid = (tl + tr) / 2;

        int left = first_active(tree[v].l, ql, qr, tl, mid);
        if(left != -1)return left;

        return first_active(tree[v].r, ql, qr, mid + 1, tr);
    }

    int last_active(int v, int ql, int qr, int tl, int tr) {
        if(!v || tree[v].sum <= 0 || qr < tl || tr < ql)return -1;
        if(tl == tr)return tl;

        int mid = (tl + tr) / 2;

        int right = last_active(tree[v].r, ql, qr, mid + 1, tr);
        if(right != -1)return right;

        return last_active(tree[v].l, ql, qr, tl, mid);
    }

    int merge_trees(int a, int b, int tl, int tr) {
        if(!a)return b;
        if(!b)return a;

        int cur = timer++;

        if(tl == tr) {
            tree[cur] = {tree[a].sum + tree[b].sum, 0, 0};
            return cur;
        }

        int mid = (tl + tr) / 2;

        int left = merge_trees(tree[a].l, tree[b].l, tl, mid);
        int right = merge_trees(tree[a].r, tree[b].r, mid + 1, tr);

        tree[cur] = join(left, right);

        return cur;
    }

public:
    PersistentSegtree() {
        n = 0;
        timer = 1;
    }

    PersistentSegtree(int _n, int mx_nodes) {
        init(_n, mx_nodes);
    }

    void init(int _n, int mx_nodes) {
        n = _n;
        timer = 1;
        tree.assign(mx_nodes, {});
    }

    // builds empty root with all zeros
    int build_empty() {
        return build_empty(0, n - 1);
    }

    // builds root from initial array, 0-indexed, size n
    int build(const vector<ll>&arr) {
        return build(0, n - 1, arr);
    }

    // creates new version: a[pos]+=val
    int add(int root, int pos,ll val) {
        return add(root, pos, val, 0, n - 1);
    }

    // creates new version: a[pos]=val
    int set(int root, int pos,ll val) {
        return set_value(root, pos, val, 0, n - 1);
    }

    // sum/count in one version over [l,r]
    ll query(int root, int l, int r) {
        if(l > r)return 0;
        l = max(l, 0);
        r = min(r, n - 1);
        if(l > r)return 0;
        return query(root, l, r, 0, n - 1);
    }

    // query difference of two prefix versions over [l,r]
    ll query(int rroot, int lroot, int l, int r) {
        if(l > r)return 0;
        return query(rroot, l, r) - query(lroot, l, r);
    }

    // value at position pos in root
    ll point_value(int root, int pos) {
        return query(root, pos, pos);
    }

    // kth active/frequency position in one root
    int kth(int root,ll k) {
        return kth_single(root, k, 0, n - 1);
    }

    // kth active/frequency position in rroot-lroot
    int kth(int rroot, int lroot,ll k) {
        return kth_range(rroot, lroot, k, 0, n - 1);
    }

    // kth active/frequency position on tree path
    int kth(int uroot, int vroot, int lcroot, int plcroot,ll k) {
        return kth_path(uroot, vroot, lcroot, plcroot, k, 0, n - 1);
    }

    // query value-id interval [l,r] on tree path
    ll path_query(int uroot, int vroot, int lcroot, int plcroot, int l, int r) {
        if(l > r)return 0;
        l = max(l, 0);
        r = min(r, n - 1);
        if(l > r)return 0;
        return path_query(uroot, vroot, lcroot, plcroot, l, r, 0, n - 1);
    }

    // first position with positive contribution in [l,r]
    int first_active(int root, int l, int r) {
        if(l > r)return -1;
        return first_active(root, l, r, 0, n - 1);
    }

    // last position with positive contribution in [l,r]
    int last_active(int root, int l, int r) {
        if(l > r)return -1;
        return last_active(root, l, r, 0, n - 1);
    }

    // merges two roots and returns new root
    int merge(int root1, int root2) {
        return merge_trees(root1, root2, 0, n - 1);
    }

    // copies root to a new root id
    int add_copy(int root) {
        int cur = timer++;
        tree[cur] = tree[root];
        return cur;
    }

    // total sum/frequency in root
    ll total(int root) {
        return tree[root].sum;
    }

    // used nodes count
    int used_nodes() {
        return timer;
    }
};

//====================================================
// Compression Helpers - 0-indexed compressed ids
//====================================================

vector<ll> comp;

void compress_values(vector<ll> v) {
    comp = v;
    sort(all(comp));
    comp.erase(unique(all(comp)), comp.end());
}

int get_id(ll x) {
    return lower_bound(all(comp), x) - comp.begin();
}

int get_leq_id(ll x) {
    return upper_bound(all(comp), x) - comp.begin() - 1;
}

int get_less_id(ll x) {
    return lower_bound(all(comp), x) - comp.begin() - 1;
}

int get_geq_id(ll x) {
    return lower_bound(all(comp), x) - comp.begin();
}

int get_greater_id(ll x) {
    return upper_bound(all(comp), x) - comp.begin();
}

bool exists_value(ll x) {
    int id = get_id(x);
    return id >= 0 && id < (int) comp.size() && comp[id] == x;
}

//====================================================
// Array Prefix Version Helpers
//====================================================

// builds roots[i+1] from roots[i] by adding compressed a[i]
vector<int> build_prefix_frequency_roots(vector<ll>&a, PersistentSegtree&pst) {
    int n = a.size();

    vector<int> roots(n + 1);
    roots[0] = pst.build_empty();

    for(int i = 0;i < n;i++) {
        int id = get_id(a[i]);
        roots[i + 1] = pst.add(roots[i], id, 1);
    }

    return roots;
}

// frequency of exact value x in a[l..r], 0-indexed l,r
ll count_equal(vector<int>&roots, PersistentSegtree&pst, int l, int r,ll x) {
    if(!exists_value(x))return 0;

    int id = get_id(x);

    return pst.query(roots[r + 1], roots[l], id, id);
}

// count values <= x in a[l..r], 0-indexed l,r
ll count_leq(vector<int>&roots, PersistentSegtree&pst, int l, int r,ll x) {
    int id = get_leq_id(x);

    if(id < 0)return 0;

    return pst.query(roots[r + 1], roots[l], 0, id);
}

// count values < x in a[l..r], 0-indexed l,r
ll count_less(vector<int>&roots, PersistentSegtree&pst, int l, int r,ll x) {
    int id = get_less_id(x);

    if(id < 0)return 0;

    return pst.query(roots[r + 1], roots[l], 0, id);
}

// count values >= x in a[l..r], 0-indexed l,r
ll count_geq(vector<int>&roots, PersistentSegtree&pst, int l, int r,ll x) {
    int id = get_geq_id(x);

    if(id >= (int) comp.size())return 0;

    return pst.query(roots[r + 1], roots[l], id, (int) comp.size() - 1);
}

// count values > x in a[l..r], 0-indexed l,r
ll count_greater(vector<int>&roots, PersistentSegtree&pst, int l, int r,ll x) {
    int id = get_greater_id(x);

    if(id >= (int) comp.size())return 0;

    return pst.query(roots[r + 1], roots[l], id, (int) comp.size() - 1);
}

// count values in [x,y] inside a[l..r], 0-indexed l,r
ll count_between(vector<int>&roots, PersistentSegtree&pst, int l, int r,ll x,ll y) {
    if(x > y)return 0;

    int left = get_geq_id(x);
    int right = get_leq_id(y);

    if(left > right)return 0;
    if(left >= (int) comp.size() || right < 0)return 0;

    left = max(left, 0);
    right = min(right, (int) comp.size() - 1);

    return pst.query(roots[r + 1], roots[l], left, right);
}

// kth smallest value in a[l..r], 0-indexed l,r
ll kth_smallest(vector<int>&roots, PersistentSegtree&pst, int l, int r, int k) {
    int id = pst.kth(roots[r + 1], roots[l], k);
    return comp[id];
}

// kth largest value in a[l..r], 0-indexed l,r
ll kth_largest(vector<int>&roots, PersistentSegtree&pst, int l, int r, int k) {
    int len = r - l + 1;
    return kth_smallest(roots, pst, l, r, len - k + 1);
}

// lower median value in a[l..r], 0-indexed l,r
ll range_median(vector<int>&roots, PersistentSegtree&pst, int l, int r) {
    int len = r - l + 1;
    return kth_smallest(roots, pst, l, r, (len + 1) / 2);
}

// minimum value in a[l..r], 0-indexed l,r
ll range_min_value(vector<int>&roots, PersistentSegtree&pst, int l, int r) {
    return kth_smallest(roots, pst, l, r, 1);
}

// maximum value in a[l..r], 0-indexed l,r
ll range_max_value(vector<int>&roots, PersistentSegtree&pst, int l, int r) {
    return kth_smallest(roots, pst, l, r, r - l + 1);
}

//====================================================
// Static Distinct Helpers
//====================================================

// prev occurrence method
// stores prev[i] as 0..n, where 0 means no previous occurrence
// root[i+1] contains prev values for indices [0..i]
vector<int> build_distinct_prev_roots(vector<ll>&a, PersistentSegtree&pst) {
    int n = a.size();

    vector<int> roots(n + 1);
    roots[0] = pst.build_empty();

    map<ll, int> last;

    for(int i = 0;i < n;i++) {
        int prev = last.count(a[i]) ? last[a[i]] + 1 : 0;
        roots[i + 1] = pst.add(roots[i], prev, 1);
        last[a[i]] = i;
    }

    return roots;
}

// distinct count in a[l..r], 0-indexed l,r
ll distinct_count_prev(vector<int>&roots, PersistentSegtree&pst, int l, int r) {
    return pst.query(roots[r + 1], roots[l], 0, l);
}

// checks if all elements in a[l..r] are distinct, static
bool all_distinct_static(vector<int>&roots, PersistentSegtree&pst, int l, int r) {
    return distinct_count_prev(roots, pst, l, r) == r - l + 1;
}

// number of duplicate positions in a[l..r]
ll duplicate_positions_count(vector<int>&roots, PersistentSegtree&pst, int l, int r) {
    return (r - l + 1) - distinct_count_prev(roots, pst, l, r);
}

// latest active occurrence method
// root[i+1] keeps only latest occurrence active in prefix [0..i]
vector<int> build_latest_active_roots(vector<ll>&a, PersistentSegtree&pst) {
    int n = a.size();

    vector<int> roots(n + 1);
    roots[0] = pst.build_empty();

    map<ll, int> last;

    for(int i = 0;i < n;i++) {
        roots[i + 1] = roots[i];

        if(last.count(a[i])) {
            roots[i + 1] = pst.add(roots[i + 1], last[a[i]], -1);
        }

        roots[i + 1] = pst.add(roots[i + 1], i, 1);
        last[a[i]] = i;
    }

    return roots;
}

// distinct count using latest active method
ll distinct_count_latest(vector<int>&roots, PersistentSegtree&pst, int l, int r) {
    return pst.query(roots[r + 1], l, r);
}

//====================================================
// Versioned Array Helpers
//====================================================

// builds initial version root from array
int build_versioned_array(PersistentSegtree&pst, vector<ll>&a) {
    return pst.build(a);
}

// creates new version by setting a[pos]=val
int new_set_version(PersistentSegtree&pst, int root, int pos,ll val) {
    return pst.set(root, pos, val);
}

// creates new version by adding delta to a[pos]
int new_add_version(PersistentSegtree&pst, int root, int pos,ll delta) {
    return pst.add(root, pos, delta);
}

// range sum in version root
ll version_range_sum(PersistentSegtree&pst, int root, int l, int r) {
    return pst.query(root, l, r);
}

// point value in version root
ll version_point_value(PersistentSegtree&pst, int root, int pos) {
    return pst.point_value(root, pos);
}

//====================================================
// Inversion Count Helper
//====================================================

// counts inversions in static array
ll count_inversions(vector<ll> a) {
    compress_values(a);

    PersistentSegtree pst((int) comp.size(), ((int) a.size() + 5) * 25);
    vector<int> roots = build_prefix_frequency_roots(a, pst);

    ll inv = 0;

    for(int i = 0;i < (int) a.size();i++) {
        int id = get_id(a[i]);

        if(id + 1 < (int) comp.size()) {
            inv += pst.query(roots[i], id + 1, (int) comp.size() - 1);
        }
    }

    return inv;
}

//====================================================
// Tree LCA Helper for Path Persistent Segment Tree
//====================================================

struct LCA {
    int n, lg, timer;
    vector<int> dep, tin, tout;
    vector<vector<int> > up;
    vector<vector<int> > g;

    LCA() {
    }

    LCA(int _n) {
        init(_n);
    }

    void init(int _n) {
        n = _n;
        lg = 1;
        while((1 << lg) <= n)lg++;

        timer = 0;

        dep.assign(n + 1, 0);
        tin.assign(n + 1, 0);
        tout.assign(n + 1, 0);
        up.assign(n + 1, vector<int>(lg, 0));
        g.assign(n + 1, {});
    }

    void add_edge(int u, int v) {
        g[u].pb(v);
        g[v].pb(u);
    }

    void dfs(int u, int p) {
        tin[u] = ++timer;
        up[u][0] = p;
        dep[u] = dep[p] + 1;

        for(int i = 1;i < lg;i++) {
            up[u][i] = up[up[u][i - 1]][i - 1];
        }

        for(auto v:g[u]) {
            if(v == p)continue;
            dfs(v, u);
        }

        tout[u] = ++timer;
    }

    void build(int root = 1) {
        dfs(root, 0);
    }

    int jump(int u, int k) {
        for(int i = 0;i < lg;i++) {
            if(k & (1 << i))u = up[u][i];
        }

        return u;
    }

    int lca(int u, int v) {
        if(dep[u] < dep[v])swap(u, v);

        u = jump(u, dep[u] - dep[v]);

        if(u == v)return u;

        for(int i = lg - 1;i >= 0;i--) {
            if(up[u][i] != up[v][i]) {
                u = up[u][i];
                v = up[v][i];
            }
        }

        return up[u][0];
    }

    int parent(int u) {
        return up[u][0];
    }

    int dist(int u, int v) {
        int lc = lca(u, v);
        return dep[u] + dep[v] - 2 * dep[lc];
    }

    int path_vertices(int u, int v) {
        return dist(u, v) + 1;
    }
};

// builds tree path roots:
// root[u] = root[parent[u]] + value[u]
void build_tree_persistent_roots(
    int u,
    int p,
    vector<vector<int> >&g,
    vector<ll>&val,
    vector<int>&roots,
    PersistentSegtree&pst
) {
    int id = get_id(val[u]);
    roots[u] = pst.add(roots[p], id, 1);

    for(auto v:g[u]) {
        if(v == p)continue;
        build_tree_persistent_roots(v, u, g, val, roots, pst);
    }
}

// kth smallest on path u-v
ll kth_smallest_on_path(
    int u,
    int v,
    int k,
    vector<int>&roots,
    LCA&lca,
    PersistentSegtree&pst
) {
    int lc = lca.lca(u, v);
    int plc = lca.parent(lc);
    int id = pst.kth(roots[u], roots[v], roots[lc], roots[plc], k);

    return comp[id];
}

// count values <= x on path u-v
ll count_leq_on_path(
    int u,
    int v,
    ll x,
    vector<int>&roots,
    LCA&lca,
    PersistentSegtree&pst
) {
    int id = get_leq_id(x);

    if(id < 0)return 0;

    int lc = lca.lca(u, v);
    int plc = lca.parent(lc);

    return pst.path_query(roots[u], roots[v], roots[lc], roots[plc], 0, id);
}

// count values < x on path u-v
ll count_less_on_path(
    int u,
    int v,
    ll x,
    vector<int>&roots,
    LCA&lca,
    PersistentSegtree&pst
) {
    int id = get_less_id(x);

    if(id < 0)return 0;

    int lc = lca.lca(u, v);
    int plc = lca.parent(lc);

    return pst.path_query(roots[u], roots[v], roots[lc], roots[plc], 0, id);
}

// count values >= x on path u-v
ll count_geq_on_path(
    int u,
    int v,
    ll x,
    vector<int>&roots,
    LCA&lca,
    PersistentSegtree&pst
) {
    int id = get_geq_id(x);

    if(id >= (int) comp.size())return 0;

    int lc = lca.lca(u, v);
    int plc = lca.parent(lc);

    return pst.path_query(roots[u], roots[v], roots[lc], roots[plc], id, (int) comp.size() - 1);
}

// count values > x on path u-v
ll count_greater_on_path(
    int u,
    int v,
    ll x,
    vector<int>&roots,
    LCA&lca,
    PersistentSegtree&pst
) {
    int id = get_greater_id(x);

    if(id >= (int) comp.size())return 0;

    int lc = lca.lca(u, v);
    int plc = lca.parent(lc);

    return pst.path_query(roots[u], roots[v], roots[lc], roots[plc], id, (int) comp.size() - 1);
}

// count exact value x on path u-v
ll count_equal_on_path(
    int u,
    int v,
    ll x,
    vector<int>&roots,
    LCA&lca,
    PersistentSegtree&pst
) {
    if(!exists_value(x))return 0;

    int id = get_id(x);
    int lc = lca.lca(u, v);
    int plc = lca.parent(lc);

    return pst.path_query(roots[u], roots[v], roots[lc], roots[plc], id, id);
}

//====================================================
// Example solve() with modes
//====================================================

void solve() {
    int mode;
    cin >> mode;

    // mode 1: kth smallest in subarray
    if(mode == 1) {
        int n, q;
        cin >> n >> q;

        vector<ll> a(n);

        for(int i = 0;i < n;i++) {
            cin >> a[i];
        }

        compress_values(a);

        PersistentSegtree pst((int) comp.size(), (n + 5) * 25);
        vector<int> roots = build_prefix_frequency_roots(a, pst);

        while(q--) {
            int l, r, k;
            cin >> l >> r >> k;
            l--, r--;

            cout << kth_smallest(roots, pst, l, r, k) << '\n';
        }
    }

    // mode 2: count <= x in subarray
    else if(mode == 2) {
        int n, q;
        cin >> n >> q;

        vector<ll> a(n);

        for(int i = 0;i < n;i++) {
            cin >> a[i];
        }

        compress_values(a);

        PersistentSegtree pst((int) comp.size(), (n + 5) * 25);
        vector<int> roots = build_prefix_frequency_roots(a, pst);

        while(q--) {
            int l, r;
            ll x;
            cin >> l >> r >> x;
            l--, r--;

            cout << count_leq(roots, pst, l, r, x) << '\n';
        }
    }

    // mode 3: count distinct values in subarray
    else if(mode == 3) {
        int n, q;
        cin >> n >> q;

        vector<ll> a(n);

        for(int i = 0;i < n;i++) {
            cin >> a[i];
        }

        PersistentSegtree pst(n + 1, (n + 5) * 25);
        vector<int> roots = build_distinct_prev_roots(a, pst);

        while(q--) {
            int l, r;
            cin >> l >> r;
            l--, r--;

            cout << distinct_count_prev(roots, pst, l, r) << '\n';
        }
    }

    // mode 4: static all-distinct check in subarray
    else if(mode == 4) {
        int n, q;
        cin >> n >> q;

        vector<ll> a(n);

        for(int i = 0;i < n;i++) {
            cin >> a[i];
        }

        PersistentSegtree pst(n + 1, (n + 5) * 25);
        vector<int> roots = build_distinct_prev_roots(a, pst);

        while(q--) {
            int l, r;
            cin >> l >> r;
            l--, r--;

            cout << (all_distinct_static(roots, pst, l, r) ? "YES" : "NO") << '\n';
        }
    }

    // mode 5: versioned array with point set and range sum
    else if(mode == 5) {
        int n, q;
        cin >> n >> q;

        vector<ll> a(n);

        for(int i = 0;i < n;i++) {
            cin >> a[i];
        }

        PersistentSegtree pst(n, (n + q + 5) * 25);
        vector<int> roots(q + 5);

        roots[0] = pst.build(a);

        int versions = 0;

        while(q--) {
            int type;
            cin >> type;

            if(type == 1) {
                int ver, pos;
                ll val;
                cin >> ver >> pos >> val;
                pos--;

                versions++;
                roots[versions] = pst.set(roots[ver], pos, val);
            } else {
                int ver, l, r;
                cin >> ver >> l >> r;
                l--, r--;

                cout << pst.query(roots[ver], l, r) << '\n';
            }
        }
    }

    // mode 6: kth smallest on tree path
    else if(mode == 6) {
        int n, q;
        cin >> n >> q;

        vector<ll> val(n + 1), vals;

        for(int i = 1;i <= n;i++) {
            cin >> val[i];
            vals.pb(val[i]);
        }

        compress_values(vals);

        PersistentSegtree pst((int) comp.size(), (n + 5) * 25);

        LCA lca(n);

        vector<vector<int> > g(n + 1);

        for(int i = 1;i < n;i++) {
            int u, v;
            cin >> u >> v;

            g[u].pb(v);
            g[v].pb(u);

            lca.add_edge(u, v);
        }

        lca.build(1);

        vector<int> roots(n + 1, 0);
        roots[0] = pst.build_empty();

        build_tree_persistent_roots(1, 0, g, val, roots, pst);

        while(q--) {
            int u, v, k;
            cin >> u >> v >> k;

            cout << kth_smallest_on_path(u, v, k, roots, lca, pst) << '\n';
        }
    }

    // mode 7: inversion count
    else if(mode == 7) {
        int n;
        cin >> n;

        vector<ll> a(n);

        for(int i = 0;i < n;i++) {
            cin >> a[i];
        }

        cout << count_inversions(a) << '\n';
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    solve();
    return 0;
}
