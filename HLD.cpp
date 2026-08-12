```// HLD, 1-indexed (adj size n+1, node 0 unused)
//   adj: vector<vector<int>> or vector<vector<pair<int,int>>> ({to, weight})
//   EDGES = 0 -> values on nodes, EDGES = 1 -> values on edges
struct HLD {
    int n, timer = 0, EDGES;
    vector<vector<pair<int, int>>> adj;
    vector<int> in, sz, big, head, par, lvl, w; // w[v] = weight of edge (par[v], v)

    template <class T>
    HLD(const vector<vector<T>> &g, int root = 1, int EDGES = 0)
        : n((int)g.size() - 1), EDGES(EDGES), adj(g.size()),
          in(n + 1), sz(n + 1), big(n + 1), head(n + 1), par(n + 1), lvl(n + 1), w(n + 1) {
        for (int u = 1; u <= n; ++u)
            for (auto e : g[u])
                if constexpr (is_same_v<T, int>) adj[u].push_back({e, 0});
                else adj[u].push_back(e);
        lvl[root] = 1, par[root] = head[root] = root;
        dfs0(root, -1), flatin(root, -1);
    }

    void dfs0(int u, int p) {
        sz[u] = 1;
        for (auto [v, c] : adj[u]) {
            if (v == p) continue;
            par[v] = u, w[v] = c, lvl[v] = lvl[u] + 1;
            dfs0(v, u);
            sz[u] += sz[v];
            if (big[u] == 0 || sz[v] > sz[big[u]]) big[u] = v;
        }
    }

    void flatin(int u, int p) {
        in[u] = timer++;
        if (big[u]) head[big[u]] = head[u], flatin(big[u], u);
        for (auto [v, c] : adj[u]) {
            if (v == p || v == big[u]) continue;
            head[v] = v, flatin(v, u);
        }
    }

    // inclusive ranges [a, b] covering the u-v path (empty in edge mode if u == v)
    vector<pair<int, int>> path(int u, int v) {
        vector<pair<int, int>> res;
        while (true) {
            if (head[u] == head[v]) {
                if (lvl[u] < lvl[v]) swap(u, v);
                if (in[v] + EDGES <= in[u]) res.push_back({in[v] + EDGES, in[u]});
                return res;
            }
            if (lvl[head[v]] > lvl[head[u]]) swap(u, v);
            res.push_back({in[head[u]], in[u]});
            u = par[head[u]];
        }
    }

    pair<int, int> subtree(int v) { return {in[v] + EDGES, in[v] + sz[v] - 1}; }
};```