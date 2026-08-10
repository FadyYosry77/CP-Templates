
/*
    Usage:

        int n;
        cin >> n;

        SmallToLarge dsu(n);       // or: Sack dsu(n);

        vector<ll> a(n);
        for(ll&x:a)
            cin >> x;

        dsu.set_values(a);

        for(int i = 0;i < n - 1;i++) {
            int u,v;
            cin >> u >> v;
            u--,v--;
            dsu.add_edge(u,v);
        }

        dsu.build(0, 3);           // fixed_k = 3

        // dsu.distinct[u]
        // dsu.mx[u]
        // dsu.mode_sum[u]
        // dsu.exact[u]
        // dsu.pairs[u]
        // dsu.frequency_classes[u]

    Extra SmallToLarge calls:

        // vector<ll> ans = dsu.subtree_mex();
        // ans = dsu.minimum_subtree_pair_difference();
        // ans = dsu.cross_child_sum_pairs(target);
        // ans = dsu.exact_distance_pairs(k);
        // ans = dsu.node_path_xor_pairs(target);
        // ans = dsu.palindrome_mask_path_pairs(26);

    Offline Sack queries must be attached before build():

        Sack sack(n);
        // set values and add edges...

        int q0 = sack.add_value_frequency_query(u, value);
        int q1 = sack.add_exact_frequency_query(u, k);
        int q2 = sack.add_at_least_frequency_query(u, k);
        int q3 = sack.add_frequency_range_query(u, l, r);

        sack.build();

        // sack.query_answer[q0], ..., sack.query_answer[q3]

    For a custom Small-to-Large problem:
        1. Add metadata fields to SmallToLarge::State.
        2. Update them in add_block().
        3. Save the answer at the end of dfs().

    For a custom Sack problem:
        1. Add global-state fields.
        2. Update them in add_node(u, delta).
        3. Save the answer in save_answer(u).

    Recursive DFS may overflow on a chain near 2e5 nodes.
*/


#include <bits/stdc++.h>
using namespace std;

using ll = long long;

/*
    Compact Small-to-Large / Sack template.

    SmallToLarge uses one std::map state per active subtree.
    Complexity: O(n log^2 n), memory O(n) after consumed maps are cleared.

    Sack uses one global compressed frequency state.
    Complexity: O(n log n), memory O(n+q), assuming O(1) add/remove.
    Attached raw-value queries add O(q log n); frequency-threshold queries are O(q).

    Both structs calculate these answers for every subtree:
        distinct[u]  = number of distinct values
        mx[u]        = maximum value frequency
        mode_sum[u]  = sum of values having frequency mx[u]
        exact[u]     = number of values occurring exactly fixed_k times
        pairs[u]     = number of unordered equal-value node pairs
        frequency_classes[u] = number of different positive frequencies present

    SmallToLarge utilities:
        cross-child equal/sum/XOR pairs
        subtree MEX, ordered successor, minimum pair difference
        fixed/variable depth counts and depth-filtered distinct values
        exact-distance pairs grouped by LCA
        inclusive node-path XOR and palindrome-mask pairs grouped by LCA

    Sack supports attached per-subtree queries for:
        frequency of a value
        number of values occurring exactly / at least k times
        number of values whose frequencies lie in [l,r]
*/

struct SmallToLarge {
    struct State {
        map<ll, ll> cnt;
        map<ll, ll> freq_cnt;
        ll mx = 0;
        ll mode_sum = 0;
        ll exact = 0;
        ll pairs = 0;
    };

    int n = 0, root = 0;
    ll fixed_k = 1;
    vector<vector<int>> g;
    vector<ll> a;
    vector<State> st;

    vector<ll> distinct;
    vector<ll> mx;
    vector<ll> mode_sum;
    vector<ll> exact;
    vector<ll> pairs;
    vector<ll> frequency_classes;

    // Create an empty tree with n nodes. O(n)
    SmallToLarge(int _n = 0) {
        init(_n);
    }

    // Reset to an empty n-node tree. O(n)
    void init(int _n) {
        n = _n;
        root = 0;
        fixed_k = 1;
        g.assign(n, {});
        a.assign(n, 0);
        st.clear();
        distinct.clear();
        mx.clear();
        mode_sum.clear();
        exact.clear();
        pairs.clear();
        frequency_classes.clear();
    }

    // Add an undirected edge. O(1)
    void add_edge(int u, int v) {
        g[u].push_back(v);
        g[v].push_back(u);
    }

    // Set all node values. O(n)
    void set_values(const vector<ll>&_a) {
        assert((int) _a.size() == n);
        a = _a;
    }

    // Reset a consumed state. O(size of its remaining map)
    void clear(State&s) {
        s.cnt.clear();
        s.freq_cnt.clear();
        s.mx = s.mode_sum = s.exact = s.pairs = 0;
    }

    // Add c occurrences of x and update all metadata. O(log size)
    // This is the main function to modify for problem-specific metadata.
    void add_block(State&s, ll x, ll c) {
        ll old = s.cnt[x];

        if(old > 0) {
            s.freq_cnt[old]--;

            if(s.freq_cnt[old] == 0)
                s.freq_cnt.erase(old);
        }

        if(old == fixed_k)
            s.exact--;

        s.pairs += old * c + c * (c - 1) / 2;

        ll now = (s.cnt[x] += c);
        s.freq_cnt[now]++;

        if(now == fixed_k)
            s.exact++;

        if(now > s.mx) {
            s.mx = now;
            s.mode_sum = x;
        } else if(now == s.mx) {
            s.mode_sum += x;
        }
    }

    // Merge b into a and consume b. O(|small| log |large|)
    // Always swap the whole State so metadata follows its container.
    void merge(State&x, State&y) {
        if(x.cnt.size() < y.cnt.size())
            swap(x, y);

        for(auto [value, count]:y.cnt)
            add_block(x, value, count);

        clear(y);
    }

    // Build every standard subtree answer. O(n log^2 n)
    void dfs(int u, int p) {
        add_block(st[u], a[u], 1);

        for(int v:g[u]) {
            if(v == p)
                continue;

            dfs(v, u);
            merge(st[u], st[v]);
        }

        distinct[u] = st[u].cnt.size();
        mx[u] = st[u].mx;
        mode_sum[u] = st[u].mode_sum;
        exact[u] = st[u].exact;
        pairs[u] = st[u].pairs;
        frequency_classes[u] = st[u].freq_cnt.size();
    }

    // Rebuild all standard answers for one fixed k. O(n log^2 n)
    void build(int _root = 0, ll _fixed_k = 1) {
        assert(!n || (0 <= _root && _root < n));
        assert(_fixed_k >= 1);
        root = _root;
        fixed_k = _fixed_k;

        st.assign(n, {});
        distinct.assign(n, 0);
        mx.assign(n, 0);
        mode_sum.assign(n, 0);
        exact.assign(n, 0);
        pairs.assign(n, 0);
        frequency_classes.assign(n, 0);

        if(n)
            dfs(root, -1);
    }

    // Generic symmetric query-before-insert merge. O(n log^2 n)
    // need(x) returns the value that can pair with x.
    template<class F>
    void cross_dfs(int u, int p, F&need, vector<map<ll,ll>>&mp,
                   vector<ll>&ans) const {
        ll ways = 0;

        for(int v:g[u]) {
            if(v == p)
                continue;

            cross_dfs(v, u, need, mp, ans);

            if(mp[u].size() < mp[v].size())
                swap(mp[u], mp[v]);

            for(auto [x, c]:mp[v]) {
                auto it = mp[u].find(need(x));

                if(it != mp[u].end())
                    ways += c * it->second;
            }

            for(auto [x, c]:mp[v])
                mp[u][x] += c;

            mp[v].clear();
        }

        ans[u] = ways;
        mp[u][a[u]]++;
    }

    // Count equal-value pairs from different child subtrees. O(n log^2 n)
    vector<ll> cross_child_equal_pairs(int _root = 0) const {
        auto need = [](ll x) {
            return x;
        };

        vector<map<ll,ll>> mp(n);
        vector<ll> ans(n);

        if(n)
            cross_dfs(_root, -1, need, mp, ans);

        return ans;
    }

    // Count pairs from different child subtrees with sum target. O(n log^2 n)
    vector<ll> cross_child_sum_pairs(ll target, int _root = 0) const {
        auto need = [target](ll x) {
            return target - x;
        };

        vector<map<ll,ll>> mp(n);
        vector<ll> ans(n);

        if(n)
            cross_dfs(_root, -1, need, mp, ans);

        return ans;
    }

    // Count pairs from different child subtrees with XOR target. O(n log^2 n)
    vector<ll> cross_child_xor_pairs(ll target, int _root = 0) const {
        auto need = [target](ll x) {
            return target ^ x;
        };

        vector<map<ll,ll>> mp(n);
        vector<ll> ans(n);

        if(n)
            cross_dfs(_root, -1, need, mp, ans);

        return ans;
    }

    // Helper for subtree_mex(). O(n log^2 n)
    void mex_dfs(int u, int p, vector<set<ll>>&bag,
                 vector<ll>&cur, vector<ll>&ans) const {
        if(a[u] >= 0)
            bag[u].insert(a[u]);

        for(int v:g[u]) {
            if(v == p)
                continue;

            mex_dfs(v, u, bag, cur, ans);

            if(bag[u].size() < bag[v].size()) {
                swap(bag[u], bag[v]);
                swap(cur[u], cur[v]);
            }

            for(ll x:bag[v])
                bag[u].insert(x);

            bag[v].clear();
        }

        while(bag[u].count(cur[u]))
            cur[u]++;

        ans[u] = cur[u];
    }

    // MEX of nonnegative values in every subtree. O(n log^2 n)
    vector<ll> subtree_mex(int _root = 0) const {
        vector<set<ll>> bag(n);
        vector<ll> cur(n), ans(n);

        if(n)
            mex_dfs(_root, -1, bag, cur, ans);

        return ans;
    }

    // Helper for nodes_at_distance_below(). O(n log n)
    // bag[u].back() stores distance 0 from u.
    void depth_dfs(int u, int p, ll k, vector<vector<ll>>&bag,
                   vector<ll>&ans) const {
        bag[u] = {1};

        for(int v:g[u]) {
            if(v == p)
                continue;

            depth_dfs(v, u, k, bag, ans);
            bag[v].push_back(0);

            if(bag[u].size() < bag[v].size())
                swap(bag[u], bag[v]);

            for(int d = 0;d < (int) bag[v].size();d++) {
                int x = (int) bag[u].size() - 1 - d;
                int y = (int) bag[v].size() - 1 - d;
                bag[u][x] += bag[v][y];
            }

            bag[v].clear();
        }

        if(k < (ll) bag[u].size())
            ans[u] = bag[u][(int) bag[u].size() - 1 - k];
    }

    // Count nodes exactly k edges below every node. O(n log n)
    vector<ll> nodes_at_distance_below(ll k, int _root = 0) const {
        assert(k >= 0);
        vector<vector<ll>> bag(n);
        vector<ll> ans(n);

        if(n)
            depth_dfs(_root, -1, k, bag, ans);

        return ans;
    }

    // Helper for variable nodes_at_distance_below(). O(n log n)
    void variable_depth_dfs(int u, int p, const vector<int>&k,
                            vector<vector<ll>>&bag, vector<ll>&ans) const {
        bag[u] = {1};

        for(int v:g[u]) {
            if(v == p)
                continue;

            variable_depth_dfs(v, u, k, bag, ans);
            bag[v].push_back(0);

            if(bag[u].size() < bag[v].size())
                swap(bag[u], bag[v]);

            for(int d = 0;d < (int) bag[v].size();d++) {
                int x = (int) bag[u].size() - 1 - d;
                int y = (int) bag[v].size() - 1 - d;
                bag[u][x] += bag[v][y];
            }

            bag[v].clear();
        }

        if(k[u] >= 0 && k[u] < (int) bag[u].size())
            ans[u] = bag[u][(int) bag[u].size() - 1 - k[u]];
    }

    // Count nodes exactly k[u] edges below each node. O(n log n)
    vector<ll> nodes_at_distance_below(const vector<int>&k, int _root = 0) const {
        assert((int) k.size() == n);
        vector<vector<ll>> bag(n);
        vector<ll> ans(n);

        if(n)
            variable_depth_dfs(_root, -1, k, bag, ans);

        return ans;
    }

    struct GapState {
        map<ll,ll> cnt;
        ll best = (1LL << 62);
    };

    // Add c copies of x while maintaining the minimum pair difference. O(log size)
    void gap_add(GapState&s, ll x, ll c) const {
        auto it = s.cnt.lower_bound(x);

        if(it != s.cnt.end() && it->first == x) {
            s.best = 0;
            it->second += c;
            return;
        }

        if(c >= 2)
            s.best = 0;

        if(it != s.cnt.end())
            s.best = min(s.best, it->first - x);

        if(it != s.cnt.begin()) {
            auto prv = prev(it);
            s.best = min(s.best, x - prv->first);
        }

        s.cnt[x] = c;
    }

    // Helper for minimum_subtree_pair_difference(). O(n log^2 n)
    void gap_dfs(int u, int p, vector<GapState>&bag, vector<ll>&ans) const {
        gap_add(bag[u], a[u], 1);

        for(int v:g[u]) {
            if(v == p)
                continue;

            gap_dfs(v, u, bag, ans);

            if(bag[u].cnt.size() < bag[v].cnt.size())
                swap(bag[u], bag[v]);

            for(auto [x,c]:bag[v].cnt)
                gap_add(bag[u], x, c);

            bag[v].cnt.clear();
        }

        ans[u] = (bag[u].best == (1LL << 62) ? -1 : bag[u].best);
    }

    // Minimum |a[x]-a[y]| over distinct nodes in every subtree, or -1. O(n log^2 n)
    vector<ll> minimum_subtree_pair_difference(int _root = 0) const {
        vector<GapState> bag(n);
        vector<ll> ans(n, -1);

        if(n)
            gap_dfs(_root, -1, bag, ans);

        return ans;
    }

    // Helper for subtree_successor(). O(n log^2 n)
    void successor_dfs(int u, int p, const vector<ll>&need,
                       vector<set<ll>>&bag, vector<ll>&ans) const {
        bag[u].insert(a[u]);

        for(int v:g[u]) {
            if(v == p)
                continue;

            successor_dfs(v, u, need, bag, ans);

            if(bag[u].size() < bag[v].size())
                swap(bag[u], bag[v]);

            for(ll x:bag[v])
                bag[u].insert(x);

            bag[v].clear();
        }

        auto it = bag[u].lower_bound(need[u]);
        ans[u] = (it == bag[u].end() ? -1 : *it);
    }

    // Smallest subtree value >= need[u], or -1. O(n log^2 n)
    vector<ll> subtree_successor(const vector<ll>&need, int _root = 0) const {
        assert((int) need.size() == n);
        vector<set<ll>> bag(n);
        vector<ll> ans(n, -1);

        if(n)
            successor_dfs(_root, -1, need, bag, ans);

        return ans;
    }

    struct DepthValueState {
        map<pair<int,ll>,ll> cnt;
        map<int,ll> distinct;
    };

    // Helper for distinct_values_at_distance_below(). O(n log^2 n)
    void depth_value_dfs(int u, int p, int dep, const vector<int>&k,
                         vector<DepthValueState>&bag, vector<ll>&ans) const {
        bag[u].cnt[{dep, a[u]}] = 1;
        bag[u].distinct[dep] = 1;

        for(int v:g[u]) {
            if(v == p)
                continue;

            depth_value_dfs(v, u, dep + 1, k, bag, ans);

            if(bag[u].cnt.size() < bag[v].cnt.size())
                swap(bag[u], bag[v]);

            for(auto [key,c]:bag[v].cnt) {
                ll&now = bag[u].cnt[key];

                if(now == 0)
                    bag[u].distinct[key.first]++;

                now += c;
            }

            bag[v].cnt.clear();
            bag[v].distinct.clear();
        }

        int wanted = dep + k[u];
        auto it = bag[u].distinct.find(wanted);
        ans[u] = (k[u] < 0 || it == bag[u].distinct.end() ? 0 : it->second);
    }

    // Distinct values exactly k[u] edges below each node. O(n log^2 n)
    vector<ll> distinct_values_at_distance_below(const vector<int>&k,
                                                  int _root = 0) const {
        assert((int) k.size() == n);
        vector<DepthValueState> bag(n);
        vector<ll> ans(n);

        if(n)
            depth_value_dfs(_root, -1, 0, k, bag, ans);

        return ans;
    }

    // Helper for exact_distance_pairs(). O(n log^2 n)
    void distance_pair_dfs(int u, int p, ll k, vector<ll>&dep,
                           vector<map<ll,ll>>&mp, vector<ll>&ans) const {
        mp[u][dep[u]] = 1;
        ll ways = 0;

        for(int v:g[u]) {
            if(v == p)
                continue;

            dep[v] = dep[u] + 1;
            distance_pair_dfs(v, u, k, dep, mp, ans);

            if(mp[u].size() < mp[v].size())
                swap(mp[u], mp[v]);

            for(auto [d,c]:mp[v]) {
                ll need = k + 2 * dep[u] - d;
                auto it = mp[u].find(need);

                if(it != mp[u].end())
                    ways += c * it->second;
            }

            for(auto [d,c]:mp[v])
                mp[u][d] += c;

            mp[v].clear();
        }

        ans[u] = ways;
    }

    // Pairs with LCA u and tree distance exactly k, for every u. O(n log^2 n)
    vector<ll> exact_distance_pairs(ll k, int _root = 0) const {
        assert(k >= 0);
        vector<ll> dep(n), ans(n);
        vector<map<ll,ll>> mp(n);

        if(n)
            distance_pair_dfs(_root, -1, k, dep, mp, ans);

        return ans;
    }

    // Helper for node_path_xor_pairs(). O(n log^2 n)
    void path_xor_dfs(int u, int p, ll target, vector<ll>&pref,
                      vector<map<ll,ll>>&mp, vector<ll>&ans) const {
        pref[u] = a[u] ^ (p == -1 ? 0 : pref[p]);
        mp[u][pref[u]] = 1;
        ll ways = 0;

        for(int v:g[u]) {
            if(v == p)
                continue;

            path_xor_dfs(v, u, target, pref, mp, ans);

            if(mp[u].size() < mp[v].size())
                swap(mp[u], mp[v]);

            for(auto [x,c]:mp[v]) {
                auto it = mp[u].find(x ^ target ^ a[u]);

                if(it != mp[u].end())
                    ways += c * it->second;
            }

            for(auto [x,c]:mp[v])
                mp[u][x] += c;

            mp[v].clear();
        }

        ans[u] = ways;
    }

    // Pairs with LCA u whose inclusive node-value path XOR is target. O(n log^2 n)
    vector<ll> node_path_xor_pairs(ll target, int _root = 0) const {
        vector<ll> pref(n), ans(n);
        vector<map<ll,ll>> mp(n);

        if(n)
            path_xor_dfs(_root, -1, target, pref, mp, ans);

        return ans;
    }

    // Helper for palindrome_mask_path_pairs(). O(bits*n log^2 n)
    void palindrome_dfs(int u, int p, int bits, vector<ll>&pref,
                        vector<map<ll,ll>>&mp, vector<ll>&ans) const {
        pref[u] = a[u] ^ (p == -1 ? 0 : pref[p]);
        mp[u][pref[u]] = 1;
        ll ways = 0;

        for(int v:g[u]) {
            if(v == p)
                continue;

            palindrome_dfs(v, u, bits, pref, mp, ans);

            if(mp[u].size() < mp[v].size())
                swap(mp[u], mp[v]);

            for(auto [x,c]:mp[v]) {
                ll base = x ^ a[u];
                auto it = mp[u].find(base);

                if(it != mp[u].end())
                    ways += c * it->second;

                for(int bit = 0;bit < bits;bit++) {
                    it = mp[u].find(base ^ (1LL << bit));

                    if(it != mp[u].end())
                        ways += c * it->second;
                }
            }

            for(auto [x,c]:mp[v])
                mp[u][x] += c;

            mp[v].clear();
        }

        ans[u] = ways;
    }

    // Pairs with LCA u whose node masks XOR to <=1 set bit. O(bits*n log^2 n)
    // Use a[u] = 1LL<<(letter-'a') for palindrome-permutable path strings.
    vector<ll> palindrome_mask_path_pairs(int bits = 26, int _root = 0) const {
        assert(0 <= bits && bits <= 62);
        vector<ll> pref(n), ans(n);
        vector<map<ll,ll>> mp(n);

        if(n)
            palindrome_dfs(_root, -1, bits, pref, mp, ans);

        return ans;
    }
};

struct Sack {
    enum QueryType {
        VALUE_FREQUENCY,
        EXACT_FREQUENCY,
        AT_LEAST_FREQUENCY,
        FREQUENCY_RANGE
    };

    struct Query {
        QueryType type;
        ll x, y;
        int id;
    };

    int n = 0, root = 0, timer = 0;
    ll fixed_k = 1;

    vector<vector<int>> g;
    vector<ll> a;
    vector<ll> values;
    vector<int> id;

    vector<int> sz, heavy, tin, tout, euler;
    vector<ll> cnt, freq_cnt, at_least_cnt, sum_by_freq;

    vector<vector<Query>> queries;
    vector<ll> query_answer;

    ll cur_distinct = 0;
    ll cur_mx = 0;
    ll cur_pairs = 0;
    ll cur_frequency_classes = 0;

    vector<ll> distinct;
    vector<ll> mx;
    vector<ll> mode_sum;
    vector<ll> exact;
    vector<ll> pairs;
    vector<ll> frequency_classes;

    // Create an empty tree with n nodes. O(n)
    Sack(int _n = 0) {
        init(_n);
    }

    // Reset to an empty n-node tree. O(n)
    void init(int _n) {
        n = _n;
        root = timer = 0;
        fixed_k = 1;
        g.assign(n, {});
        a.assign(n, 0);
        values.clear();
        id.clear();
        sz.clear();
        heavy.clear();
        tin.clear();
        tout.clear();
        euler.clear();
        cnt.clear();
        freq_cnt.clear();
        at_least_cnt.clear();
        sum_by_freq.clear();
        queries.assign(n, {});
        query_answer.clear();
        distinct.clear();
        mx.clear();
        mode_sum.clear();
        exact.clear();
        pairs.clear();
        frequency_classes.clear();
    }

    // Add an undirected edge. O(1)
    void add_edge(int u, int v) {
        g[u].push_back(v);
        g[v].push_back(u);
    }

    // Set all node values. O(n)
    void set_values(const vector<ll>&_a) {
        assert((int) _a.size() == n);
        a = _a;
    }

    // Attach query: frequency of value x in subtree u. Returns query ID. O(1)
    int add_value_frequency_query(int u, ll x) {
        assert(0 <= u && u < n);
        int qid = query_answer.size();
        queries[u].push_back({VALUE_FREQUENCY, x, 0, qid});
        query_answer.push_back(0);
        return qid;
    }

    // Attach query: number of values occurring exactly k times. Returns query ID. O(1)
    int add_exact_frequency_query(int u, ll k) {
        assert(0 <= u && u < n);
        int qid = query_answer.size();
        queries[u].push_back({EXACT_FREQUENCY, k, 0, qid});
        query_answer.push_back(0);
        return qid;
    }

    // Attach query: number of values occurring at least k times. Returns query ID. O(1)
    int add_at_least_frequency_query(int u, ll k) {
        assert(0 <= u && u < n);
        int qid = query_answer.size();
        queries[u].push_back({AT_LEAST_FREQUENCY, k, 0, qid});
        query_answer.push_back(0);
        return qid;
    }

    // Attach query: number of values with frequency in [l,r]. Returns query ID. O(1)
    int add_frequency_range_query(int u, ll l, ll r) {
        assert(0 <= u && u < n);
        int qid = query_answer.size();
        queries[u].push_back({FREQUENCY_RANGE, l, r, qid});
        query_answer.push_back(0);
        return qid;
    }

    // Compress values and allocate dense frequency arrays. O(n log n)
    void compress() {
        values = a;
        sort(values.begin(), values.end());
        values.erase(unique(values.begin(), values.end()), values.end());

        id.resize(n);

        for(int i = 0;i < n;i++)
            id[i] = lower_bound(values.begin(), values.end(), a[i]) - values.begin();

        cnt.assign(values.size(), 0);
        freq_cnt.assign(n + 1, 0);
        at_least_cnt.assign(n + 2, 0);
        sum_by_freq.assign(n + 1, 0);
    }

    // Build subtree sizes, Euler order, and heavy children. O(n)
    void dfs_size(int u, int p) {
        sz[u] = 1;
        heavy[u] = -1;
        tin[u] = timer;
        euler[timer++] = u;

        for(int v:g[u]) {
            if(v == p)
                continue;

            dfs_size(v, u);
            sz[u] += sz[v];

            if(heavy[u] == -1 || sz[v] > sz[heavy[u]])
                heavy[u] = v;
        }

        tout[u] = timer;
    }

    // Add or remove one node from the global state. Amortized O(1)
    // This is the main function to modify for problem-specific state.
    void add_node(int u, int delta) {
        assert(delta == 1 || delta == -1);
        int x = id[u];
        ll value = values[x];
        ll old = cnt[x];

        if(old > 0) {
            if(freq_cnt[old] == 1)
                cur_frequency_classes--;

            freq_cnt[old]--;
            sum_by_freq[old] -= value;
        }

        ll now = old + delta;

        if(delta == 1)
            at_least_cnt[now]++;
        else
            at_least_cnt[old]--;

        cur_pairs += now * (now - 1) / 2 - old * (old - 1) / 2;
        cnt[x] = now;

        if(!old && now)
            cur_distinct++;
        if(old && !now)
            cur_distinct--;

        if(now > 0) {
            if(freq_cnt[now] == 0)
                cur_frequency_classes++;

            freq_cnt[now]++;
            sum_by_freq[now] += value;
            cur_mx = max(cur_mx, now);
        }

        while(cur_mx > 0 && freq_cnt[cur_mx] == 0)
            cur_mx--;
    }

    // Add/remove one Euler subtree interval. O(subtree size)
    void add_subtree(int u, int delta) {
        for(int t = tin[u];t < tout[u];t++)
            add_node(euler[t], delta);
    }

    // Frequency of one raw value in the current Sack state. O(log n)
    ll current_value_frequency(ll value) const {
        auto it = lower_bound(values.begin(), values.end(), value);

        if(it == values.end() || *it != value)
            return 0;

        return cnt[it - values.begin()];
    }

    // Number of values occurring exactly k times in the current state. O(1)
    ll current_exact_frequency(ll k) const {
        if(k <= 0 || k > n)
            return 0;

        return freq_cnt[k];
    }

    // Number of values occurring at least k times in the current state. O(1)
    ll current_at_least_frequency(ll k) const {
        if(k <= 1)
            return cur_distinct;
        if(k > n)
            return 0;

        return at_least_cnt[k];
    }

    // Number of values whose current frequency is in [l,r]. O(1)
    ll current_frequency_range(ll l, ll r) const {
        l = max(1LL, l);
        r = min<ll>(n, r);

        if(l > r)
            return 0;

        return at_least_cnt[l] - at_least_cnt[r + 1];
    }

    // Save answers for u. O(1 + attached queries, with O(log n) raw-value lookup)
    // Add problem-specific answers here.
    void save_answer(int u) {
        distinct[u] = cur_distinct;
        mx[u] = cur_mx;
        mode_sum[u] = (cur_mx ? sum_by_freq[cur_mx] : 0);
        exact[u] = (fixed_k <= n ? freq_cnt[fixed_k] : 0);
        pairs[u] = cur_pairs;
        frequency_classes[u] = cur_frequency_classes;

        for(const Query&q:queries[u]) {
            if(q.type == VALUE_FREQUENCY)
                query_answer[q.id] = current_value_frequency(q.x);
            else if(q.type == EXACT_FREQUENCY)
                query_answer[q.id] = current_exact_frequency(q.x);
            else if(q.type == AT_LEAST_FREQUENCY)
                query_answer[q.id] = current_at_least_frequency(q.x);
            else
                query_answer[q.id] = current_frequency_range(q.x, q.y);
        }
    }

    // Sack DFS. O(n log n) total node additions/removals.
    void dfs(int u, int p, bool keep) {
        for(int v:g[u]) {
            if(v != p && v != heavy[u])
                dfs(v, u, false);
        }

        if(heavy[u] != -1)
            dfs(heavy[u], u, true);

        for(int v:g[u]) {
            if(v != p && v != heavy[u])
                add_subtree(v, 1);
        }

        add_node(u, 1);
        save_answer(u);

        if(!keep)
            add_subtree(u, -1);
    }

    // Build answers. O(n log n + q log n) with q attached raw-value queries.
    void build(int _root = 0, ll _fixed_k = 1) {
        assert(!n || (0 <= _root && _root < n));
        assert(_fixed_k >= 1);
        root = _root;
        fixed_k = _fixed_k;
        timer = 0;

        compress();

        sz.assign(n, 0);
        heavy.assign(n, -1);
        tin.assign(n, 0);
        tout.assign(n, 0);
        euler.assign(n, 0);

        distinct.assign(n, 0);
        mx.assign(n, 0);
        mode_sum.assign(n, 0);
        exact.assign(n, 0);
        pairs.assign(n, 0);
        frequency_classes.assign(n, 0);

        cur_distinct = cur_mx = cur_pairs = cur_frequency_classes = 0;
        fill(query_answer.begin(), query_answer.end(), 0);

        if(n) {
            dfs_size(root, -1);
            dfs(root, -1, true);
        }
    }
};
