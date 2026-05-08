#define potassium() (ios_base::sync_with_stdio(false), cin.tie(NULL));

#include <bits/stdc++.h>
// #include <ext/pb_ds/assoc_container.hpp>
// #include <ext/pb_ds/tree_policy.hpp>
const char nl = '\n';
#define ll long long
#define int long long
#define tc                                                                     \
  ll _;                                                                        \
  cin >> _;                                                                    \
  while (_--)
// #define int ll;
using namespace std;

// template<class T>
// using Tree = tree<T, null_type, less<T>, rb_tree_tag, tree_order_statistics_node_update>;
void TXT(string name = "") {
#ifndef ONLINE_JUDGE
    if (name.empty()) {
        freopen("input.txt", "r", stdin);
        freopen("output.txt", "w", stdout);
    } else {
        freopen((name + ".in").c_str(), "r", stdin);
        freopen((name + ".out").c_str(), "w", stdout);
    }
#else
    if (!name.empty()) {
      freopen((name + ".in").c_str(), "r", stdin);
      freopen((name + ".out").c_str(), "w", stdout);
    }
#endif
}
// #pragma GCC optimize("O3,unroll-loops,fast-math")
// #pragma GCC target("avx2,sse4.2,fma")
typedef vector<ll> vll;
typedef vector<vll> mat;
typedef map<ll, ll> mll;
typedef pair<ll, ll> pll;
typedef double ld;
typedef unsigned long long ull;
#define all(x) x.begin(), x.end()
#define rall(x) (x).rbegin(), (x).rend()
// #define sz(x) ((long long)x.size())
#define F first
#define S second
#define pb push_back
#define mem(a,v) memset(a, v, sizeof(a))
#define fr front()
#define bk back()
#define pr make_pair
#define eb emplace_back
#define dbg(x) cerr << #x << " = " << (x) << nl
template<class T>
void _print(const vector<T> &v) {
    //.\go.bat
    cerr << "[ ";
    for (auto &x : v) { _print(x); cerr << " "; }
    cerr << "]";
}
using vi = vector<int>;
#define rep(i,a,b) for (int i=(a); i<(b); ++i)
int dx[4]={2,-2,0,0};
int dy[4]={0,0,2,-2};
const ld eps = 1e-9;
const ld PI = acos(-1);
mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
const ll inf = 4e18, mod = 1e9 + 7; // 998244353
struct LCA {
    int n, lg, timer;
    vector<int> dep, sz, tin, tout;
    vector<vector<int> > par;
    vector<vector<ll> > mx, mn, sum;
    vector<vector<pair<int, ll> > > g;

    LCA() {
    }

    LCA(int _n) {
        init(_n);
    }

    // initializes all arrays
    void init(int _n) {
        n = _n;
        lg = 1;
        while((1 << lg) <= n)lg++;

        timer = 0;
        dep.assign(n + 1, 0);
        sz.assign(n + 1, 0);
        tin.assign(n + 1, 0);
        tout.assign(n + 1, 0);

        par.assign(n + 1, vector<int>(lg, 0));
        mx.assign(n + 1, vector<ll>(lg, 0));
        mn.assign(n + 1, vector<ll>(lg, 4e18));
        sum.assign(n + 1, vector<ll>(lg, 0));

        g.assign(n + 1, {});
    }

    // adds an undirected weighted edge
    void add_edge(int u, int v, ll w = 1) {
        g[u].pb({v, w});
        g[v].pb({u, w});
    }

    // clears the tree edges only
    void clear_edges() {
        for(int i = 1;i <= n;i++)g[i].clear();
    }

    // dfs preprocessing for binary lifting
    void dfs(int u, int p, ll w) {
        tin[u] = ++timer;
        par[u][0] = p;
        dep[u] = dep[p] + 1;
        sz[u] = 1;

        mx[u][0] = w;
        mn[u][0] = w;
        sum[u][0] = w;

        for(int i = 1;i < lg;i++) {
            par[u][i] = par[par[u][i - 1]][i - 1];
            mx[u][i] = max(mx[u][i - 1], mx[par[u][i - 1]][i - 1]);
            mn[u][i] = min(mn[u][i - 1], mn[par[u][i - 1]][i - 1]);
            sum[u][i] = sum[u][i - 1] + sum[par[u][i - 1]][i - 1];
        }

        for(auto [v,c]:g[u]) {
            if(v == p)continue;
            dfs(v, u, c);
            sz[u] += sz[v];
        }

        tout[u] = ++timer;
    }

    // builds the LCA structure from a root
    void build(int root = 1) {
        dfs(root, 0, 0);
        mn[root][0] = 4e18;
    }

    // returns true if u is ancestor of v
    bool is_ancestor(int u, int v) {
        return tin[u] <= tin[v] && tout[v] <= tout[u];
    }

    // returns true if u is inside subtree of v
    bool in_subtree(int u, int v) {
        return is_ancestor(v, u);
    }

    // moves u up by k edges
    int jump(int u, int k) {
        for(int i = 0;i < lg;i++) {
            if(k & (1 << i))u = par[u][i];
        }
        return u;
    }

    // returns the kth ancestor of u
    int kth_ancestor(int u, int k) {
        return jump(u, k);
    }

    // returns lca of u and v
    int lca(int u, int v) {
        if(dep[u] < dep[v])swap(u, v);

        u = jump(u, dep[u] - dep[v]);

        if(u == v)return u;

        for(int i = lg - 1;i >= 0;i--) {
            if(par[u][i] != par[v][i]) {
                u = par[u][i];
                v = par[v][i];
            }
        }

        return par[u][0];
    }

    // returns distance in edges between u and v
    int dist(int u, int v) {
        int l = lca(u, v);
        return dep[u] + dep[v] - 2 * dep[l];
    }

    // returns number of vertices on path u-v
    int path_vertices_count(int u, int v) {
        return dist(u, v) + 1;
    }

    // returns number of edges on path u-v
    int path_edges_count(int u, int v) {
        return dist(u, v);
    }

    // returns kth node from u to v, 0th node is u
    int go(int u, int v, int k) {
        int l = lca(u, v);
        int d1 = dep[u] - dep[l];
        int d2 = dep[v] - dep[l];
        int d = d1 + d2;

        assert(k <= d);

        if(k <= d1)return jump(u, k);

        k -= d1;
        return jump(v, d2 - k);
    }

    // returns true if x lies on path u-v
    bool on_path(int x, int u, int v) {
        return dist(u, x) + dist(x, v) == dist(u, v);
    }

    // returns true if paths a-b and c-d intersect
    bool paths_intersect(int a, int b, int c, int d) {
        int x = lca(a, b);
        int y = lca(c, d);
        return on_path(x, c, d) || on_path(y, a, b);
    }

    // returns the child of ancestor anc that contains node u
    int child_under(int anc, int u) {
        assert(is_ancestor(anc, u) && anc != u);
        return jump(u, dep[u] - dep[anc] - 1);
    }

    // returns subtree size of u with original root
    int subtree_size(int u) {
        return sz[u];
    }

    // returns subtree size of u if the tree is rooted at root
    int subtree_size_rooted(int u, int root) {
        if(u == root)return n;

        if(!is_ancestor(u, root))return sz[u];

        int c = child_under(u, root);
        return n - sz[c];
    }

    // returns parent of u
    int parent(int u) {
        return par[u][0];
    }

    // returns depth of u
    int depth(int u) {
        return dep[u] - 1;
    }

    // compares nodes by dfs order
    bool by_tin(int u, int v) {
        return tin[u] < tin[v];
    }

    // sorts nodes by dfs order
    void sort_by_tin(vector<int>&v) {
        sort(all(v), [&](int a, int b) {
            return tin[a] < tin[b];
        });
    }

    // removes duplicated nodes from a vector
    void unique_nodes(vector<int>&v) {
        sort(all(v));
        v.erase(unique(all(v)), v.end());
    }

    // returns lca when the tree is considered rooted at r
    int reroot_lca(int u, int v, int r) {
        int a = lca(u, v);
        int b = lca(u, r);
        int c = lca(v, r);

        if(a == b)return c;
        if(a == c)return b;
        return a;
    }

    // returns maximum edge weight on path u-v
    ll get_max_edge(int u, int v) {
        ll ans = 0;

        if(dep[u] < dep[v])swap(u, v);

        int diff = dep[u] - dep[v];

        for(int i = lg - 1;i >= 0;i--) {
            if(diff & (1 << i)) {
                ans = max(ans, mx[u][i]);
                u = par[u][i];
            }
        }

        if(u == v)return ans;

        for(int i = lg - 1;i >= 0;i--) {
            if(par[u][i] != par[v][i]) {
                ans = max(ans, mx[u][i]);
                ans = max(ans, mx[v][i]);
                u = par[u][i];
                v = par[v][i];
            }
        }

        ans = max(ans, mx[u][0]);
        ans = max(ans, mx[v][0]);

        return ans;
    }

    // returns minimum edge weight on path u-v
    ll get_min_edge(int u, int v) {
        ll ans = 4e18;

        if(dep[u] < dep[v])swap(u, v);

        int diff = dep[u] - dep[v];

        for(int i = lg - 1;i >= 0;i--) {
            if(diff & (1 << i)) {
                ans = min(ans, mn[u][i]);
                u = par[u][i];
            }
        }

        if(u == v)return ans;

        for(int i = lg - 1;i >= 0;i--) {
            if(par[u][i] != par[v][i]) {
                ans = min(ans, mn[u][i]);
                ans = min(ans, mn[v][i]);
                u = par[u][i];
                v = par[v][i];
            }
        }

        ans = min(ans, mn[u][0]);
        ans = min(ans, mn[v][0]);

        return ans;
    }

    // returns sum of edge weights on path u-v
    ll get_sum_edge(int u, int v) {
        ll ans = 0;

        if(dep[u] < dep[v])swap(u, v);

        int diff = dep[u] - dep[v];

        for(int i = lg - 1;i >= 0;i--) {
            if(diff & (1 << i)) {
                ans += sum[u][i];
                u = par[u][i];
            }
        }

        if(u == v)return ans;

        for(int i = lg - 1;i >= 0;i--) {
            if(par[u][i] != par[v][i]) {
                ans += sum[u][i];
                ans += sum[v][i];
                u = par[u][i];
                v = par[v][i];
            }
        }

        ans += sum[u][0];
        ans += sum[v][0];

        return ans;
    }

    // returns max, min, and sum of edges on path u-v
    array<ll, 3> get_path_info(int u, int v) {
        ll a = get_max_edge(u, v);
        ll b = get_min_edge(u, v);
        ll c = get_sum_edge(u, v);
        return {a, b, c};
    }

    // returns farthest node from src among given nodes
    int farthest_node(int src, vector<int> nodes) {
        int best = src;

        for(auto u:nodes) {
            if(dist(src, u) > dist(src, best))best = u;
        }

        return best;
    }

    // returns diameter endpoints among given nodes
    pair<int, int> diameter(vector<int> nodes) {
        int a = nodes[0];
        int b = farthest_node(a, nodes);
        int c = farthest_node(b, nodes);
        return {b, c};
    }

    // returns diameter length among given nodes
    int diameter_length(vector<int> nodes) {
        auto [a,b] = diameter(nodes);
        return dist(a, b);
    }

    // returns compressed list containing nodes and their LCAs for virtual tree
    vector<int> virtual_tree_nodes(vector<int> nodes) {
        sort_by_tin(nodes);

        int k = nodes.size();

        for(int i = 0;i + 1 < k;i++) {
            nodes.pb(lca(nodes[i], nodes[i + 1]));
        }

        sort_by_tin(nodes);
        nodes.erase(unique(all(nodes)), nodes.end());

        return nodes;
    }
};
void fady(){
    ll n;cin>>n;
    LCA lca(n);
    for(int i=0;i<n-1;i++){
        ll u,v,w;cin>>u>>v>>w;
        lca.add_edge(u,v,w);
    }
    lca.build(1);
    ll q;cin>>q;
    while(q--){
        ll u,v;cin>>u>>v;
        cout<<"u = "<<u<<", v = "<<v<<nl;
        cout<<"lca = "<<lca.lca(u,v)<<nl;
        cout<<"distance = "<<lca.dist(u,v)<<nl;
        cout<<"max edge = "<<lca.get_max_edge(u,v)<<nl;
        cout<<"min edge = "<<lca.get_min_edge(u,v)<<nl;
        cout<<"sum edge = "<<lca.get_sum_edge(u,v)<<nl;
        cout<<"middle node from u to v = "<<lca.go(u,v,lca.dist(u,v)/2)<<nl;
        cout<<nl;
    }
}
signed main() {
    TXT();
    potassium();
    fady();
    // just normal comment for nonsense readers...
    return 0;
}
