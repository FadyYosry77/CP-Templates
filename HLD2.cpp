// HEAVY-LIGHT DECOMPOSITION -- C++17
// N = number of vertices; every path has O(log N) heavy-chain pieces.
//
// # Build                                        O(N)
// # Space                                        O(N)
// # dfs0 / flatin                                O(N) total during build
// # path                                         O(log N) ranges
// # path_ordered                                 O(log N) directed ranges
// # for_each_path / for_each_path_ordered        O(log N) callbacks
// # fold_path                                    O(log N * range-query cost)
// # fold_path_ordered                            O(log N * range-query cost)
// # subtree                                      O(1)
// # subtree_size / position / vertex_at          O(1)
// # is_ancestor                                  O(1)
// # lca                                          O(log N)
// # lca_rooted / is_on_path                      O(log N)
// # distance / weighted_distance                 O(log N)
// # path_length                                  O(log N)
// # kth_ancestor                                 O(log N)
// # kth_on_path / jump                           O(log N)
// # child_on_path                                O(log N)
// # rooted_subtree                               O(log N), at most two ranges
// # edge_node / edge_position                    O(1)
// # flatten_nodes / flatten_edges                O(N)
// # restore_nodes                                O(N)
//
// HOW TO RECOGNIZE AN HLD PROBLEM
//
// 1. The input is a static tree and queries/updates concern an arbitrary path
//    u..v: path sum/min/max/gcd/xor, recoloring, adding, assigning, or finding
//    the first position satisfying a condition.
//
// 2. You already know how to solve the operation on one ARRAY interval with a
//    segment tree/Fenwick/sparse table. HLD turns a tree path into O(log N)
//    contiguous intervals in one flattened array.
//
// 3. The operation is on a subtree as well as paths. HLD's DFS order makes
//    every rooted subtree one contiguous interval.
//
// 4. Values live on EDGES. Store each edge at the position of its deeper
//    endpoint and construct HLD with EDGES=1; the LCA position is then skipped.
//
// 5. Direction matters: strings, matrices, function composition, hashes, or
//    searching from u toward v. Use path_ordered/fold_path_ordered and respect
//    each part's reversed flag.
//
// WHEN NOT TO USE IT
//
// - Only LCA/distances/kth ancestor are needed: binary lifting may be simpler.
// - Only offline path additions/counting are needed: tree difference + DFS.
// - The tree topology changes: use a link-cut tree or another dynamic-tree tool.
// - Every query is a subtree query only: Euler tour alone is enough.
//
// CONVENTIONS
//
// - Tree is 1-indexed; adjacency size is N+1 and vertex 0 is unused.
// - Flattened ranges are 0-indexed and inclusive.
// - EDGES=0 stores node values; EDGES=1 stores edge (par[v],v) at in[v].
// - path(u,u) is one point in node mode and empty in edge mode.
// - path() may return pieces in any order: use it for commutative operations.
// - path_ordered() returns pieces in traversal order from u to v.
// - Recursive DFS is concise; for a pathological chain near stack limits,
//   increase stack size or convert the two build DFS routines to iterative DFS.

#include <bits/stdc++.h>
using namespace std;

using ll = long long;

// HLD, 1-indexed (adj size n+1, node 0 unused)
//   adj: vector<vector<int>> or vector<vector<pair<int,int>>> ({to, weight})
//   EDGES = 0 -> values on nodes, EDGES = 1 -> values on edges
struct HLD {
    // Inclusive flattened interval. If reversed=true, traverse r,r-1,...,l.
    struct Part {
        int l, r;
        bool reversed;
    };

    int n, timer = 0, EDGES, rt;
    vector<vector<pair<int, int> > > adj;
    vector<int> in, out, rev, sz, big, head, par, lvl, w;
    vector<ll> root_weight; // weighted distance from construction root

    template<class T>
    HLD(const vector<vector<T> >&g, int root = 1, int edges = 0)
        : n((int) g.size() - 1), EDGES(edges), rt(root), adj(g.size()),
          in(n + 1), out(n + 1), rev(n), sz(n + 1), big(n + 1),
          head(n + 1), par(n + 1), lvl(n + 1), w(n + 1), root_weight(n + 1) {
        assert(n>=1&&1<=root&&root<=n&&(EDGES==0||EDGES==1));

        for(int u = 1;u <= n;u++) {
            for(auto e:g[u]) {
                if constexpr(is_same_v<T, int>)adj[u].push_back({e, 0});
                else adj[u].push_back(e);
            }
        }

        lvl[root] = 1;
        par[root] = head[root] = root;
        dfs0(root, -1);
        flatin(root, -1);
    }

    // Finds subtree sizes, parents, depths, weights and heavy children. O(N)
    void dfs0(int u, int p) {
        sz[u] = 1;

        for(auto [v,c]:adj[u]) {
            if(v == p)continue;

            par[v] = u;
            w[v] = c;
            lvl[v] = lvl[u] + 1;
            root_weight[v] = root_weight[u] + c;

            dfs0(v, u);
            sz[u] += sz[v];

            if(big[u] == 0 || sz[v] > sz[big[u]])big[u] = v;
        }
    }

    // Assigns heavy-first Euler positions and chain heads. O(N)
    void flatin(int u, int p) {
        in[u] = timer;
        rev[timer] = u;
        timer++;

        if(big[u]) {
            head[big[u]] = head[u];
            flatin(big[u], u);
        }

        for(auto [v,c]:adj[u]) {
            (void) c;
            if(v == p || v == big[u])continue;
            head[v] = v;
            flatin(v, u);
        }

        out[u] = in[u] + sz[u] - 1;
    }

    // Inclusive ranges covering u-v; order is irrelevant. O(log N)
    vector<pair<int, int> > path(int u, int v) const {
        check_vertex(u);
        check_vertex(v);
        vector<pair<int, int> > answer;

        while(true) {
            if(head[u] == head[v]) {
                if(lvl[u] < lvl[v])swap(u, v);

                if(in[v] + EDGES <= in[u])
                    answer.push_back({in[v] + EDGES, in[u]});

                return answer;
            }

            if(lvl[head[v]] > lvl[head[u]])swap(u, v);
            answer.push_back({in[head[u]], in[u]});
            u = par[head[u]];
        }
    }

    // Directed ranges in exact u-to-v traversal order. O(log N)
    vector<Part> path_ordered(int u, int v) const {
        check_vertex(u);
        check_vertex(v);

        vector<Part> left, right;

        while(head[u] != head[v]) {
            if(lvl[head[u]] >= lvl[head[v]]) {
                left.push_back({in[head[u]], in[u], true});
                u = par[head[u]];
            } else {
                right.push_back({in[head[v]], in[v], false});
                v = par[head[v]];
            }
        }

        if(lvl[u] >= lvl[v]) {
            if(in[v] + EDGES <= in[u])
                left.push_back({in[v] + EDGES, in[u], true});
        } else {
            if(in[u] + EDGES <= in[v])
                left.push_back({in[u] + EDGES, in[v], false});
        }

        reverse(right.begin(), right.end());
        left.insert(left.end(), right.begin(), right.end());
        return left;
    }

    // Calls work(l,r) for every path range; commutative use. O(log N)
    template<class Work>
    void for_each_path(int u, int v, Work work) const {
        for(auto [l,r]:path(u, v))work(l, r);
    }

    // Calls work(l,r,reversed) in exact u-to-v order. O(log N)
    template<class Work>
    void for_each_path_ordered(int u, int v, Work work) const {
        for(Part part:path_ordered(u, v))
            work(part.l, part.r, part.reversed);
    }

    // Folds commutative range answers over a path. O(log N * query)
    template<class Value, class Query, class Merge>
    Value fold_path(int u, int v, Value identity, Query query, Merge merge) const {
        Value answer = identity;

        for(auto [l,r]:path(u, v))
            answer = merge(answer, query(l, r));

        return answer;
    }

    // Folds direction-sensitive range answers from u toward v.
    // query(l,r,reversed) must return the interval in requested direction.
    template<class Value, class Query, class Merge>
    Value fold_path_ordered(
        int u, int v, Value identity, Query query, Merge merge) const {
        Value answer = identity;

        for(Part part:path_ordered(u, v))
            answer = merge(answer, query(part.l, part.r, part.reversed));

        return answer;
    }

    // Flattened subtree range; may be empty for an edge-mode leaf. O(1)
    pair<int, int> subtree(int v) const {
        check_vertex(v);
        return {in[v] + EDGES, out[v]};
    }

    // Number of vertices in v's construction-root subtree. O(1)
    int subtree_size(int v) const {
        check_vertex(v);
        return sz[v];
    }

    // Flattened position of a vertex. O(1)
    int position(int v) const {
        check_vertex(v);
        return in[v];
    }

    // Vertex stored at one flattened position. O(1)
    int vertex_at(int position) const {
        assert(0<=position&&position<n);
        return rev[position];
    }

    // Whether ancestor is an ancestor of vertex under the build root. O(1)
    bool is_ancestor(int ancestor, int vertex) const {
        check_vertex(ancestor);
        check_vertex(vertex);
        return in[ancestor] <= in[vertex] && out[vertex] <= out[ancestor];
    }

    // Lowest common ancestor. O(log N)
    int lca(int u, int v) const {
        check_vertex(u);
        check_vertex(v);

        while(head[u] != head[v]) {
            if(lvl[head[u]] < lvl[head[v]])swap(u, v);
            u = par[head[u]];
        }

        return lvl[u] < lvl[v] ? u : v;
    }

    // LCA when the tree is conceptually rerooted at new_root. O(log N)
    int lca_rooted(int u, int v, int new_root) const {
        check_vertex(new_root);
        int a = lca(u, v);
        int b = lca(u, new_root);
        int c = lca(v, new_root);
        int answer = a;
        if(lvl[b] > lvl[answer])answer = b;
        if(lvl[c] > lvl[answer])answer = c;
        return answer;
    }

    // Number of edges on u-v. O(log N)
    int distance(int u, int v) const {
        int ancestor = lca(u, v);
        return lvl[u] + lvl[v] - 2 * lvl[ancestor];
    }

    // Number of represented values on u-v in the selected mode. O(log N)
    int path_length(int u, int v) const {
        return distance(u, v) + (EDGES == 0);
    }

    // Whether x lies on the simple u-v path. O(log N)
    bool is_on_path(int x, int u, int v) const {
        check_vertex(x);
        return distance(u, x) + distance(x, v) == distance(u, v);
    }

    // Sum of input edge weights on u-v. O(log N)
    ll weighted_distance(int u, int v) const {
        int ancestor = lca(u, v);
        return root_weight[u] + root_weight[v] - 2 * root_weight[ancestor];
    }

    // k-th ancestor, with k=0 returning v; 0 means nonexistent. O(log N)
    int kth_ancestor(int v, int k) const {
        check_vertex(v);
        if(k < 0 || k >= lvl[v])return 0;

        while(k > 0) {
            int chain_nodes = lvl[v] - lvl[head[v]];

            if(k <= chain_nodes)return rev[in[v] - k];

            k -= chain_nodes + 1;
            v = par[head[v]];
        }

        return v;
    }

    // k-th vertex on u-v, zero-based from u; 0 means invalid. O(log N)
    int kth_on_path(int u, int v, int k) const {
        int ancestor = lca(u, v);
        int up = lvl[u] - lvl[ancestor];
        int down = lvl[v] - lvl[ancestor];
        if(k < 0 || k > up + down)return 0;
        return k <= up ? kth_ancestor(u, k) : kth_ancestor(v, up + down - k);
    }

    // Alias for moving k edges from u toward v. O(log N)
    int jump(int u, int v, int k) const {
        return kth_on_path(u, v, k);
    }

    // Child of ancestor on the route to descendant; 0 if invalid. O(log N)
    int child_on_path(int ancestor, int descendant) const {
        check_vertex(ancestor);
        check_vertex(descendant);
        if(ancestor == descendant || !is_ancestor(ancestor, descendant))return 0;
        return kth_ancestor(descendant, lvl[descendant] - lvl[ancestor] - 1);
    }

    // Subtree of v when the conceptual root is new_root, NODE MODE only.
    // Returns 0..2 disjoint inclusive flattened ranges. O(log N)
    vector<pair<int, int> > rooted_subtree(int v, int new_root) const {
        check_vertex(v);
        check_vertex(new_root);
        assert(EDGES==0);

        if(v == new_root)return {{0, n - 1}};
        if(!is_ancestor(v, new_root))return {{in[v], out[v]}};

        int child = child_on_path(v, new_root);
        vector<pair<int, int> > answer;

        if(in[child] > 0)answer.push_back({0, in[child] - 1});
        if(out[child] + 1 < n)answer.push_back({out[child] + 1, n - 1});

        return answer;
    }

    // Deeper endpoint representing edge (u,v); 0 if not adjacent. O(1)
    int edge_node(int u, int v) const {
        check_vertex(u);
        check_vertex(v);
        if(par[u] == v)return u;
        if(par[v] == u)return v;
        return 0;
    }

    // Position representing edge (u,v), or -1 if not adjacent. O(1)
    int edge_position(int u, int v) const {
        int child = edge_node(u, v);
        return child ? in[child] : -1;
    }

    // Converts 1-indexed node values into flattened order. O(N)
    template<class Value>
    vector<Value> flatten_nodes(const vector<Value>&value) const {
        assert((int)value.size()==n+1);
        vector<Value> flat(n);

        for(int v = 1;v <= n;v++)flat[in[v]] = value[v];

        return flat;
    }

    // Converts parent-edge weights into flattened order. O(N)
    // Position in[root] receives root_value because the root has no parent edge.
    template<class Value>
    vector<Value> flatten_edges(Value root_value = Value()) const {
        vector<Value> flat(n, root_value);

        for(int v = 1;v <= n;v++) {
            if(v != rt)flat[in[v]] = (Value) w[v];
        }

        return flat;
    }

    // Converts a flattened array back to a 1-indexed node array. O(N)
    template<class Value>
    vector<Value> restore_nodes(const vector<Value>&flat) const {
        assert((int)flat.size()==n);
        vector<Value> value(n + 1);

        for(int position = 0;position < n;position++)
            value[rev[position]] = flat[position];

        return value;
    }

    // Explicit path vertices from u to v; intended for debugging/output. O(length)
    vector<int> path_vertices(int u, int v) const {
        vector<int> answer;

        for(Part part:path_ordered(u, v)) {
            if(part.reversed) {
                for(int position = part.r;position >= part.l;position--)
                    answer.push_back(rev[position]);
            } else {
                for(int position = part.l;position <= part.r;position++)
                    answer.push_back(rev[position]);
            }
        }

        return answer;
    }

    // Explicit edge child-vertices from u to v; intended for debugging. O(length)
    vector<int> path_edge_nodes(int u, int v) const {
        assert(EDGES==1);
        return path_vertices(u, v);
    }

    // Validates a public vertex argument. O(1)
    void check_vertex(int v) const {
        assert(1<=v&&v<=n);
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Tree: 1 connects to 2,3; 2 connects to 4,5; 3 connects to 6.
    vector<vector<int> > tree(7);
    auto add_edge = [&](int u, int v) {
        tree[u].push_back(v);
        tree[v].push_back(u);
    };

    add_edge(1, 2);
    add_edge(1, 3);
    add_edge(2, 4);
    add_edge(2, 5);
    add_edge(3, 6);

    HLD nodes(tree, 1, 0);
    vector<ll> value = {0, 10, 20, 30, 40, 50, 60};
    vector<ll> flat = nodes.flatten_nodes(value);

    auto range_sum = [&](int l, int r) {
        ll answer = 0;
        for(int position = l;position <= r;position++)answer += flat[position];
        return answer;
    };

    ll path_sum = nodes.fold_path(
        4, 6, 0LL, range_sum, [](ll a, ll b) { return a + b; }
    );

    assert(path_sum==40+20+10+30+60);
    assert(nodes.lca(4,5)==2);
    assert(nodes.distance(4,6)==4);
    assert(nodes.kth_on_path(4,6,2)==1);
    assert(nodes.path_vertices(4,6)==vector<int>({4,2,1,3,6}));
    assert(nodes.subtree(2)==make_pair(nodes.in[2],nodes.in[2]+2));

    string labels = " abcdef";
    string path_string = nodes.fold_path_ordered(
        4, 6, string(),
        [&](int l, int r, bool reversed) {
            string part;
            if(reversed) {
                for(int position = r;position >= l;position--)
                    part.push_back(labels[nodes.vertex_at(position)]);
            } else {
                for(int position = l;position <= r;position++)
                    part.push_back(labels[nodes.vertex_at(position)]);
            }
            return part;
        },
        [](string a, const string&b) { return a + b; }
    );

    assert(path_string=="dbacf");

    vector<vector<pair<int, int> > > weighted(7);
    auto add_weighted_edge = [&](int u, int v, int weight) {
        weighted[u].push_back({v, weight});
        weighted[v].push_back({u, weight});
    };

    add_weighted_edge(1, 2, 5);
    add_weighted_edge(1, 3, 7);
    add_weighted_edge(2, 4, 11);
    add_weighted_edge(2, 5, 13);
    add_weighted_edge(3, 6, 17);

    HLD edges(weighted, 1, 1);
    vector<int> edge_flat = edges.flatten_edges<int>();

    int path_weight = edges.fold_path(
        4, 6, 0,
        [&](int l, int r) {
            int answer = 0;
            for(int position = l;position <= r;position++)
                answer += edge_flat[position];
            return answer;
        },
        [](int a, int b) { return a + b; }
    );

    assert(path_weight==11+5+7+17);
    assert(edges.weighted_distance(4,6) == path_weight);
    assert(edges.edge_position(2,5) == edges.in[5]);
    assert(edges.path(3,3).empty());

    cout << "node path sum: " << path_sum << '\n';
    cout << "ordered path: " << path_string << '\n';
    cout << "edge path sum: " << path_weight << '\n';
}
