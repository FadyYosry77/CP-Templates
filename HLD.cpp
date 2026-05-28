struct HLD {
    int root = 0;
    vector<vector<int>> g;
    explicit HLD(int n) : g(n) { }
    void add(int u, int v) {
        g[u].push_back(v);
        g[v].push_back(u);
    }
    vector<int> &operator[](int u) {
        return g[u];
    }
    int timer = 0;
    vector<int> tin, tout, depth, sz, head, par, vertex;
    void init(int rt = 0) {
        root = rt;
        tin = tout = depth = head = par = vertex = vector<int>(g.size());
        sz.resize(g.size(), 1);
        par[root] = head[root] = root;
        dfs(root); dfs2(root);
    }
    void dfs(int u){
        for(int & v : g[u]){
            depth[v] = depth[u] + 1;
            par[v] = u;
            g[v].erase(find(g[v].begin(), g[v].end(), u));
            dfs(v);
            sz[u] += sz[v];
            if(sz[v] > sz[g[u][0]]) swap(v, g[u][0]);
        }
    }
    void dfs2(int u){
        tin[u] = timer++;
        vertex[tin[u]] = u;
        for(int v : g[u]){
            head[v] = v == g[u][0] ? head[u]: v;
            dfs2(v);
        }
        tout[u] = timer - 1;
    }
    int jump(int u, int k) {
        if(k > depth[u]) return -1;
        int d = depth[u] - k;
        while (depth[head[u]] > d) { u = par[head[u]]; }
        return vertex[tin[u] - depth[u] + d];
    }
    bool isAnc(int u, int v) {
        return tin[u] <= tin[v] && tin[v] <= tout[u];
    }
    int lca(int u, int v) {
        if(depth[u] > depth[v]) swap(u, v);
        if(isAnc(u, v)) return u;
        while(head[u] != head[v]){
            if (depth[head[u]] > depth[head[v]]) u = par[head[u]];
            else v = par[head[v]];
        }
        return depth[u] < depth[v] ? u : v;
    }
    int dis(int u, int v) {
        return depth[u] + depth[v] - 2 * depth[lca(u, v)];
    }
    auto path_sorted(int u, int v) {
        vector<tuple<int, int, bool>> U, V; // from, to, reversed
        while(head[u] != head[v]){
            if(depth[head[u]] > depth[head[v]]){
                U.emplace_back(tin[head[u]], tin[u], true);
                u = par[head[u]];
            }else{
                V.emplace_back(tin[head[v]], tin[v], false);
                v = par[head[v]];
            }
        }
        if(depth[u] < depth[v]) V.emplace_back(tin[u], tin[v], false);
        else U.emplace_back(tin[v], tin[u], true);
        U.insert(U.end(), V.rbegin(), V.rend());
        return U;
    }
};