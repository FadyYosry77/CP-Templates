struct SuffixAutomaton {
    struct State {
        array<int, 26> next;
        int link = -1, len = 0;
        bool is_clone = false;
        int first_pos = 0;
        long long cnt = 0;
        vector<int> inv_link;

        State() {
            next.fill(-1);
        }
    };

    vector<State> st;
    string s;
    int last = 0;

    // Build SAM. O(n)
    SuffixAutomaton(const string&_s = "", bool occ = true, bool inv = true) {
        build(_s, occ, inv);
    }

    // Convert char to index. O(1)
    int id(char c) const {
        return c - 'a';
    }

    // Rebuild SAM. O(n)
    void build(const string&_s, bool occ = true, bool inv = true) {
        s = _s;
        st.clear();
        st.reserve(max(1ll, 2 * (int) s.size()));
        st.emplace_back();
        last = 0;

        for(char c:s) extend(c);

        if(occ) build_occurrences();
        if(inv) build_inv_links();
    }

    // Add one character. Amortized O(1)
    void extend(char c) {
        int x = id(c);
        int cur = st.size();
        st.emplace_back();

        st[cur].len = st[last].len + 1;
        st[cur].first_pos = st[cur].len - 1;
        st[cur].cnt = 1;

        int p = last;

        while(p >= 0 && st[p].next[x] == -1) {
            st[p].next[x] = cur;
            p = st[p].link;
        }

        if(p == -1) {
            st[cur].link = 0;
        } else {
            int q = st[p].next[x];

            if(st[p].len + 1 == st[q].len) {
                st[cur].link = q;
            } else {
                int clone = st.size();
                st.push_back(st[q]);

                st[clone].len = st[p].len + 1;
                st[clone].is_clone = true;
                st[clone].cnt = 0;

                while(p >= 0 && st[p].next[x] == q) {
                    st[p].next[x] = clone;
                    p = st[p].link;
                }

                st[q].link = st[cur].link = clone;
            }
        }

        last = cur;
    }

    // States by len decreasing. O(S+n)
    vector<int> order_desc() const {
        int mx = 0,S = st.size();

        for(auto&v:st)
            mx = max(mx, v.len);

        vector<int> cnt(mx + 1), ord(S);

        for(auto&v:st)
            cnt[v.len]++;

        for(int i = 1;i <= mx;i++)
            cnt[i] += cnt[i - 1];

        for(int i = S - 1;i >= 0;i--)
            ord[--cnt[st[i].len]] = i;

        reverse(ord.begin(), ord.end());

        return ord;
    }

    // Build occurrence counts. O(S+n)
    void build_occurrences() {
        for(int i = 0;i < (int) st.size();i++)
            st[i].cnt = (i && !st[i].is_clone);

        for(int v:order_desc()) {
            if(st[v].link != -1)
                st[st[v].link].cnt += st[v].cnt;
        }
    }

    // Build suffix-link tree. O(S)
    void build_inv_links() {
        for(auto&v:st)
            v.inv_link.clear();

        for(int i = 1;i < (int) st.size();i++)
            st[st[i].link].inv_link.push_back(i);
    }

    // State after walking t. O(|t|)
    int go_state(const string&t) const {
        int v = 0;

        for(char c:t) {
            int x = id(c);

            if(x < 0 || x >= 26 || st[v].next[x] == -1)
                return -1;

            v = st[v].next[x];
        }

        return v;
    }

    // Check substring. O(|t|)
    bool contains(const string&t) const {
        return go_state(t) != -1;
    }

    // Count distinct substrings. O(S)
    long long count_distinct_substrings() const {
        long long ans = 0;

        for(int i = 1;i < (int) st.size();i++)
            ans += st[i].len - st[st[i].link].len;

        return ans;
    }

    // Total length of distinct substrings. O(S)
    long long total_length_distinct_substrings() const {
        long long ans = 0;

        for(int i = 1;i < (int) st.size();i++) {
            long long l = st[st[i].link].len + 1;
            long long r = st[i].len;
            long long c = r - l + 1;

            ans += c * (l + r) / 2;
        }

        return ans;
    }

    // Distinct substrings by length. O(S+n)
    vector<long long> distinct_by_length() const {
        int n = s.size();
        vector<long long> diff(n + 2), ans(n + 1);

        for(int i = 1;i < (int) st.size();i++) {
            int l = st[st[i].link].len + 1;
            int r = st[i].len;

            diff[l]++;
            diff[r + 1]--;
        }

        for(int i = 1;i <= n;i++) {
            diff[i] += diff[i - 1];
            ans[i] = diff[i];
        }

        return ans;
    }

    // kth lexicographical substring, 1-indexed. O(S+26*ans)
    string kth_substring(long long k) const {
        const long long INF = 4000000000000000000LL;
        int S = st.size();
        vector<long long> dp(S, -1);

        function<long long(int)> dfs = [&](int v) {
            if(dp[v] != -1)
                return dp[v];

            long long res = 0;

            for(int c = 0;c < 26;c++) {
                int u = st[v].next[c];

                if(u != -1) {
                    res += 1 + dfs(u);
                    res = min(res, INF);
                }
            }

            return dp[v] = res;
        };

        dfs(0);

        if(k <= 0 || k > dp[0])
            return "";

        string ans;
        int v = 0;

        while(k) {
            for(int c = 0;c < 26;c++) {
                int u = st[v].next[c];

                if(u == -1)
                    continue;

                long long block = 1 + dp[u];

                if(k <= block) {
                    ans.push_back(char('a' + c));
                    k--;

                    if(!k)
                        return ans;

                    v = u;
                    break;
                }

                k -= block;
            }
        }

        return ans;
    }

    // Count occurrences of t. O(|t|)
    // Requires build_occurrences().
    long long count_occurrences(const string&t) const {
        int v = go_state(t);
        return v == -1 ? 0 : st[v].cnt;
    }

    // First occurrence of t. O(|t|)
    int first_occurrence(const string&t) const {
        int v = go_state(t);

        if(v == -1)
            return -1;

        return st[v].first_pos - (int) t.size() + 1;
    }

    // Helper for all_occurrences. O(subtree)
    void report_all(int v,int len, vector<int>&res) const {
        if(!st[v].is_clone)
            res.push_back(st[v].first_pos - len + 1);

        for(int u:st[v].inv_link)
            report_all(u, len, res);
    }

    // All occurrences of t. O(|t|+subtree+occ log occ)
    // Requires build_inv_links().
    vector<int> all_occurrences(const string&t) const {
        int v = go_state(t);

        if(v == -1)
            return {};

        vector<int> res;

        report_all(v, t.size(), res);

        sort(res.begin(), res.end());

        return res;
    }

    // Longest substring appearing >= k times. O(S)
    // Requires build_occurrences().
    string longest_substring_occurring_at_least(long long k) const {
        int best = 0, pos = -1;

        for(int i = 1;i < (int) st.size();i++) {
            if(st[i].cnt >= k && st[i].len > best) {
                best = st[i].len;
                pos = st[i].first_pos;
            }
        }

        if(!best)
            return "";

        return s.substr(pos - best + 1, best);
    }

    // Longest repeated substring. O(S)
    // Requires build_occurrences().
    string longest_repeated_substring() const {
        return longest_substring_occurring_at_least(2);
    }

    // Max len * occurrences. O(S)
    // Requires build_occurrences().
    long long max_repeat_value() const {
        long long ans = 0;

        for(int i = 1;i < (int) st.size();i++)
            ans = max(ans, 1LL * st[i].len * st[i].cnt);

        return ans;
    }

    // Max occurrence count for every length. O(S log S+n log S)
    // Requires build_occurrences().
    vector<long long> max_occurrences_by_length() const {
        int n = s.size();
        vector<vector<pair<int, long long> > > add(n + 2);
        vector<long long> ans(n + 1);

        for(int i = 1;i < (int) st.size();i++) {
            int l = st[st[i].link].len + 1;
            int r = st[i].len;

            add[l].push_back({r, st[i].cnt});
        }

        priority_queue<pair<long long,int> > pq;

        for(int len = 1;len <= n;len++) {
            for(auto&[r,c]:add[len])
                pq.push({c, r});

            while(!pq.empty() && pq.top().second < len)
                pq.pop();

            ans[len] = pq.empty() ? 0 : pq.top().first;
        }

        return ans;
    }

    // Longest repeated substring without overlap. O(S+n)
    string longest_non_overlapping_repeat() const {
        int S = st.size();
        const int INF = 1e9;

        vector<int> mn(S, INF), mx(S, -INF);

        for(int i = 1;i < S;i++) {
            if(!st[i].is_clone) {
                mn[i] = st[i].first_pos;
                mx[i] = st[i].first_pos;
            }
        }

        for(int v:order_desc()) {
            if(st[v].link != -1) {
                int p = st[v].link;

                mn[p] = min(mn[p], mn[v]);
                mx[p] = max(mx[p], mx[v]);
            }
        }

        int best = 0, state = -1;

        for(int i = 1;i < S;i++) {
            int can = min(st[i].len, mx[i] - mn[i]);

            if(can > best) {
                best = can;
                state = i;
            }
        }

        if(!best)
            return "";

        return s.substr(st[state].first_pos - best + 1, best);
    }

    // Shortest absent over first k letters. O(S*k)
    string shortest_absent(int k) const {
        int S = st.size();
        vector<int> dp(S, -1);

        function<int(int)> dfs = [&](int v) {
            if(dp[v] != -1)
                return dp[v];

            int res = 1e9;

            for(int c = 0;c < k;c++) {
                int u = st[v].next[c];

                if(u == -1) {
                    res = 1;
                    break;
                }

                res = min(res, 1 + dfs(u));
            }

            return dp[v] = res;
        };

        dfs(0);

        string ans;
        int v = 0;

        while(true) {
            if(dp[v] == 1) {
                for(int c = 0;c < k;c++) {
                    if(st[v].next[c] == -1) {
                        ans.push_back(char('a' + c));
                        return ans;
                    }
                }
            }

            for(int c = 0;c < k;c++) {
                int u = st[v].next[c];

                if(u != -1 && dp[v] == 1 + dp[u]) {
                    ans.push_back(char('a' + c));
                    v = u;
                    break;
                }
            }
        }
    }

    // Minimal cyclic shift. O(n*26)
    static string minimal_cyclic_shift(const string&s) {
        SuffixAutomaton sa(s + s, false, false);

        string ans;
        int v = 0, n = s.size();

        for(int i = 0;i < n;i++) {
            for(int c = 0;c < 26;c++) {
                int u = sa.st[v].next[c];

                if(u != -1) {
                    ans.push_back(char('a' + c));
                    v = u;
                    break;
                }
            }
        }

        return ans;
    }

    // LCS of two strings. O(|A|+|B|)
    static string LCS(const string&A, const string&B) {
        SuffixAutomaton sa(A, false, false);

        int v = 0, len = 0, best = 0, bestpos = 0;

        for(int i = 0;i < (int) B.size();i++) {
            int c = B[i] - 'a';

            while(v && (c < 0 || c >= 26 || sa.st[v].next[c] == -1)) {
                v = sa.st[v].link;
                len = sa.st[v].len;
            }

            if(c >= 0 && c < 26 && sa.st[v].next[c] != -1) {
                v = sa.st[v].next[c];
                len++;
            } else {
                v = 0;
                len = 0;
            }

            if(len > best) {
                best = len;
                bestpos = i;
            }
        }

        return B.substr(bestpos - best + 1, best);
    }

    // LCS with positions in A and B. O(|A|+|B|)
    static tuple<string,int,int> LCS_with_pos(const string&A, const string&B) {
        SuffixAutomaton sa(A, false, false);

        int v = 0, len = 0, best = 0, posA = 0, posB = 0;

        for(int i = 0;i < (int) B.size();i++) {
            int c = B[i] - 'a';

            while(v && (c < 0 || c >= 26 || sa.st[v].next[c] == -1)) {
                v = sa.st[v].link;
                len = sa.st[v].len;
            }

            if(c >= 0 && c < 26 && sa.st[v].next[c] != -1) {
                v = sa.st[v].next[c];
                len++;
            } else {
                v = 0;
                len = 0;
            }

            if(len > best) {
                best = len;
                posB = i - best + 1;
                posA = sa.st[v].first_pos - best + 1;
            }
        }

        return {B.substr(posB, best), posA, posB};
    }

    // Count distinct common substrings of A and B. O(A+B+S)
    static long long count_common_distinct(const string&A, const string&B) {
        SuffixAutomaton sa(A, false, false);

        int S = sa.st.size();
        vector<int> best(S);

        int v = 0, len = 0;

        for(char ch:B) {
            int c = ch - 'a';

            while(v && (c < 0 || c >= 26 || sa.st[v].next[c] == -1)) {
                v = sa.st[v].link;
                len = sa.st[v].len;
            }

            if(c >= 0 && c < 26 && sa.st[v].next[c] != -1) {
                v = sa.st[v].next[c];
                len++;
            } else {
                v = 0;
                len = 0;
            }

            best[v] = max(best[v], len);
        }

        for(int x:sa.order_desc()) {
            if(sa.st[x].link != -1) {
                int p = sa.st[x].link;
                best[p] = max(best[p], min(best[x], sa.st[p].len));
            }
        }

        long long ans = 0;

        for(int i = 1;i < S;i++) {
            int low = sa.st[sa.st[i].link].len;
            ans += max(0ll, best[i] - low);
        }

        return ans;
    }

    // LCS among many strings. O(total+K*S)
    static string LCS_many(vector<string> ss) {
        if(ss.empty())
            return "";

        int id = 0;

        for(int i = 1;i < (int) ss.size();i++) {
            if(ss[i].size() < ss[id].size())
                id = i;
        }

        swap(ss[0], ss[id]);

        SuffixAutomaton sa(ss[0], false, false);

        int S = sa.st.size();
        vector<int> common(S);

        for(int i = 0;i < S;i++)
            common[i] = sa.st[i].len;

        for(int idx = 1;idx < (int) ss.size();idx++) {
            vector<int> cur(S);

            int v = 0, len = 0;

            for(char ch:ss[idx]) {
                int c = ch - 'a';

                while(v && (c < 0 || c >= 26 || sa.st[v].next[c] == -1)) {
                    v = sa.st[v].link;
                    len = sa.st[v].len;
                }

                if(c >= 0 && c < 26 && sa.st[v].next[c] != -1) {
                    v = sa.st[v].next[c];
                    len++;
                } else {
                    v = 0;
                    len = 0;
                }

                cur[v] = max(cur[v], len);
            }

            for(int x:sa.order_desc()) {
                if(sa.st[x].link != -1) {
                    int p = sa.st[x].link;
                    cur[p] = max(cur[p], min(cur[x], sa.st[p].len));
                }
            }

            for(int i = 0;i < S;i++)
                common[i] = min(common[i], cur[i]);
        }

        int best = 0, state = 0;

        for(int i = 1;i < S;i++) {
            if(common[i] > best) {
                best = common[i];
                state = i;
            }
        }

        if(!best)
            return "";

        return ss[0].substr(sa.st[state].first_pos - best + 1, best);
    }

    // Count substrings common to all strings. O(total+K*S)
    static long long count_common_many(vector<string> ss) {
        if(ss.empty())
            return 0;

        int id = 0;

        for(int i = 1;i < (int) ss.size();i++) {
            if(ss[i].size() < ss[id].size())
                id = i;
        }

        swap(ss[0], ss[id]);

        SuffixAutomaton sa(ss[0], false, false);

        int S = sa.st.size();
        vector<int> common(S);

        for(int i = 0;i < S;i++)
            common[i] = sa.st[i].len;

        for(int idx = 1;idx < (int) ss.size();idx++) {
            vector<int> cur(S);

            int v = 0, len = 0;

            for(char ch:ss[idx]) {
                int c = ch - 'a';

                while(v && (c < 0 || c >= 26 || sa.st[v].next[c] == -1)) {
                    v = sa.st[v].link;
                    len = sa.st[v].len;
                }

                if(c >= 0 && c < 26 && sa.st[v].next[c] != -1) {
                    v = sa.st[v].next[c];
                    len++;
                } else {
                    v = 0;
                    len = 0;
                }

                cur[v] = max(cur[v], len);
            }

            for(int x:sa.order_desc()) {
                if(sa.st[x].link != -1) {
                    int p = sa.st[x].link;
                    cur[p] = max(cur[p], min(cur[x], sa.st[p].len));
                }
            }

            for(int i = 0;i < S;i++)
                common[i] = min(common[i], cur[i]);
        }

        long long ans = 0;

        for(int i = 1;i < S;i++) {
            int low = sa.st[sa.st[i].link].len;
            ans += max(0ll, common[i] - low);
        }

        return ans;
    }
};

void solve() {
    string s;
    cin >> s;

    SuffixAutomaton sa(s);

    cout << "string: " << s << '\n';
    cout << "states: " << sa.st.size() << '\n';
    cout << "distinct substrings: " << sa.count_distinct_substrings() << '\n';
    cout << "total length of distinct substrings: " << sa.total_length_distinct_substrings() << '\n';
    cout << "longest repeated substring: " << sa.longest_repeated_substring() << '\n';
    cout << "max len * occurrences: " << sa.max_repeat_value() << '\n';
    cout << "longest non-overlapping repeat: " << sa.longest_non_overlapping_repeat() << '\n';

    cout << "\nsubstring queries:\n";

    int q;
    cin >> q;

    while(q--) {
        string t;
        cin >> t;

        cout << "pattern: " << t << '\n';
        cout << "contains: " << sa.contains(t) << '\n';
        cout << "occurrences: " << sa.count_occurrences(t) << '\n';
        cout << "first occurrence: " << sa.first_occurrence(t) << '\n';

        vector<int> occ = sa.all_occurrences(t);

        cout << "all occurrences:";
        for(auto x:occ)
            cout << " " << x;
        cout << '\n';
    }

    cout << "\nkth substring queries:\n";

    int kq;
    cin >> kq;

    while(kq--) {
        long long k;
        cin >> k;

        cout << "k = " << k << " -> " << sa.kth_substring(k) << '\n';
    }

    cout << "\nshortest absent:\n";

    int alphabet_size;
    cin >> alphabet_size;

    cout << sa.shortest_absent(alphabet_size) << '\n';

    cout << "\nminimal cyclic shift:\n";
    cout << SuffixAutomaton::minimal_cyclic_shift(s) << '\n';

    cout << "\nLCS with another string:\n";

    string b;
    cin >> b;

    cout << SuffixAutomaton::LCS(s, b) << '\n';

    auto [lcs,pos_s,pos_b] = SuffixAutomaton::LCS_with_pos(s, b);

    cout << "lcs: " << lcs << '\n';
    cout << "pos in s: " << pos_s << '\n';
    cout << "pos in b: " << pos_b << '\n';

    cout << "common distinct substrings with b: ";
    cout << SuffixAutomaton::count_common_distinct(s, b) << '\n';

    cout << "\nLCS many strings:\n";

    int m;
    cin >> m;

    vector<string> v(m);
    v[0] = s;

    for(int i = 1;i < m;i++)
        cin >> v[i];

    cout << "LCS many: " << SuffixAutomaton::LCS_many(v) << '\n';
    cout << "common substrings many: " << SuffixAutomaton::count_common_many(v) << '\n';

    cout << "\ndistinct substrings by length:\n";

    vector<long long> by_len = sa.distinct_by_length();

    for(int i = 1;i < (int) by_len.size();i++)
        cout << "len " << i << ": " << by_len[i] << '\n';

    cout << "\nmax occurrences by length:\n";

    vector<long long> mx_occ = sa.max_occurrences_by_length();

    for(int i = 1;i < (int) mx_occ.size();i++)
        cout << "len " << i << ": " << mx_occ[i] << '\n';
}
