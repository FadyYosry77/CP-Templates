// Small instructions to remember

// Think in this exact order:

// 1. What does T store?
//    Store enough information to describe one component.

// 2. What is id?
//    The value of an empty component.

// 3. merge(x,y)
//    How do two components combine at the same node?

// 4. add_edge(x)
//    What changes when that component crosses one edge?

// 5. add_node(x)
//    What does the current node itself contribute?

// 6. dp[u]
//    Final state when the whole tree is rooted at u.

template<class T, class Merge, class AddEdge, class AddVertex>
struct Reroot {
    ll n;
    vector<vector<ll> >&adj;

    T id;
    Merge merge;
    AddEdge add_edge;
    AddVertex add_vertex;

    vector<T> down;
    vector<T> up;
    vector<T> ans;

    Reroot(
        vector<vector<ll> >&adj,
        T id,
        Merge merge,
        AddEdge add_edge,
        AddVertex add_vertex
    ) : adj(adj), id(id), merge(merge), add_edge(add_edge), add_vertex(add_vertex) {
        n = adj.size();

        down.assign(n, id);
        up.assign(n, id);
        ans.assign(n, id);
    }

    void dfs1(ll u,ll p) {
        T cur = id;
        ll pe = -1;

        for(ll i = 0;i < adj[u].size();i++) {
            ll v = adj[u][i];

            if(v == p) {
                pe = i;
                continue;
            }

            dfs1(v, u);

            cur = merge(
                cur,
                add_edge(down[v], u, i)
            );
        }

        down[u] = add_vertex(cur, u, pe);
    }

    void dfs2(ll u,ll p) {
        ll m = adj[u].size();

        vector<T> a(m, id);

        for(ll i = 0;i < m;i++) {
            ll v = adj[u][i];

            if(v == p)
                a[i] = add_edge(up[u], u, i);
            else
                a[i] = add_edge(down[v], u, i);
        }

        vector<T> pre(m + 1, id);
        vector<T> suf(m + 1, id);

        for(ll i = 0;i < m;i++)
            pre[i + 1] = merge(pre[i], a[i]);

        for(ll i = m - 1;i >= 0;i--)
            suf[i] = merge(a[i], suf[i + 1]);

        ans[u] = add_vertex(pre[m], u, -1);

        for(ll i = 0;i < m;i++) {
            ll v = adj[u][i];

            if(v == p) continue;

            T without_v = merge(
                pre[i],
                suf[i + 1]
            );

            up[v] = add_vertex(without_v, u, i);

            dfs2(v, u);
        }
    }

    vector<T> run(ll root = 0) {
        dfs1(root, -1);
        dfs2(root, -1);
        return ans;
    }
};
void solve() {
    ll n; 
    cin >> n;

    vector<vector<ll>> adj(n);

    // read tree

    using T=...;

    // Empty component.
    // What should merge(x,id) return?
    T id=...;

    // Ask:
    // If two independent neighbor components
    // meet at the same node,
    // how do I combine their summaries?
    auto merge=[&](T x,T y)->T {
        return ...;
    };

    // Ask:
    // What happens to a component
    // when it crosses one edge?
    auto add_edge=[&](T x,ll u,ll idx)->T {
        return ...;
    };

    // Ask:
    // After merging all neighbor components,
    // how does node u itself affect the state?
    auto add_node=[&](T x,ll u,ll par)->T {
        return ...;
    };

    Reroot reroot(
        adj,
        id,
        merge,
        add_edge,
        add_node
    );

    vector<T> dp=reroot.run();

    // dp[u] = answer/state when u is the root

    // extract final answer
}
