// # Build Eertree                         O(n)
// # Number of nodes                       O(n)
// # Space                                 O(n)
// # add_char                              amortized O(1)
// # distinct                              O(1)
// # get_range                             O(1)
// # get_palindrome                        O(length)
// # count_occurrences with sort           O(n log n)
// # count_occurrences with reverse ids    O(n)
// # longest_palindrome                    O(length of answer)
// # all_palindromes                       O(total output size)
// # suffix_chain                          O(number of palindromic suffixes)
// # slow palindrome partition DP          O(n^2)
// # series-link partition DP              about O(n log n), fast in practice
#include <bits/stdc++.h>
using namespace std;

struct Eertree {
    static const int SIG = 26;
    static const char BASE = 'a';

    struct Node {
        int len;
        int link;
        int first_pos;

        long long raw_occ;
        long long occ;

        int diff;
        int series_link;

        int suffix_count;

        array<int, SIG> nxt;

        Node(int _len = 0) {
            len = _len;
            link = 0;
            first_pos = -1;

            raw_occ = 0;
            occ = 0;

            diff = 0;
            series_link = 0;

            suffix_count = 0;

            nxt.fill(0);
        }
    };

    vector<Node> tree;
    string s;

    int last;
    int longest_node;
    int longest_odd_node;
    int longest_even_node;

    bool occ_ready;

    vector<int> longest_suffix_node_at;
    vector<int> created_node_at;
    vector<int> distinct_after_prefix;

    Eertree() {
        reset();
    }

    // Solves: Clear the structure and start a new Eertree.
    void reset() {
        tree.clear();
        s.clear();

        tree.push_back(Node(-1)); // node 0: odd root
        tree.push_back(Node(0));  // node 1: even root

        tree[0].link = 0;
        tree[1].link = 0;

        tree[0].series_link = 0;
        tree[1].series_link = 0;

        tree[0].suffix_count = 0;
        tree[1].suffix_count = 0;

        last = 1;

        longest_node = -1;
        longest_odd_node = -1;
        longest_even_node = -1;

        occ_ready = false;

        longest_suffix_node_at.clear();
        created_node_at.clear();
        distinct_after_prefix.clear();
    }

    // Solves: Find the largest palindromic suffix that can be extended by s[pos].
    int get_link(int v, int pos) {
        while (true) {
            int left = pos - tree[v].len - 1;

            if (left >= 0 && s[left] == s[pos]) {
                return v;
            }

            v = tree[v].link;
        }
    }

    // Solves: Add one character online and return {longest suffix node, was a new palindrome created?}.
    pair<int, bool> add_char(char ch) {
        int pos = (int)s.size();
        s += ch;

        int c = ch - BASE;
        int cur = get_link(last, pos);

        bool created = false;

        if (tree[cur].nxt[c] == 0) {
            Node node(tree[cur].len + 2);
            node.first_pos = pos;

            if (node.len == 1) {
                node.link = 1;
            } else {
                int suffix_candidate = get_link(tree[cur].link, pos);
                node.link = tree[suffix_candidate].nxt[c];
            }

            node.diff = node.len - tree[node.link].len;

            if (tree[node.link].diff == node.diff) {
                node.series_link = tree[node.link].series_link;
            } else {
                node.series_link = node.link;
            }

            node.suffix_count = tree[node.link].suffix_count + 1;

            tree.push_back(node);

            int id = (int)tree.size() - 1;
            tree[cur].nxt[c] = id;

            created = true;

            if (longest_node == -1 || tree[id].len > tree[longest_node].len) {
                longest_node = id;
            }

            if (tree[id].len % 2 == 1) {
                if (longest_odd_node == -1 || tree[id].len > tree[longest_odd_node].len) {
                    longest_odd_node = id;
                }
            } else {
                if (longest_even_node == -1 || tree[id].len > tree[longest_even_node].len) {
                    longest_even_node = id;
                }
            }
        }

        last = tree[cur].nxt[c];

        tree[last].raw_occ++;

        longest_suffix_node_at.push_back(last);
        created_node_at.push_back(created ? last : -1);
        distinct_after_prefix.push_back(distinct());

        occ_ready = false;

        return {last, created};
    }

    // Solves: Build the Eertree for a full string.
    void build(const string &str) {
        for (char ch : str) {
            add_char(ch);
        }
    }

    // Solves: Return number of nodes including the two roots.
    int size() const {
        return (int)tree.size();
    }

    // Solves: Return whether the real string is empty.
    bool empty() const {
        return s.empty();
    }

    // Solves: Check if a node is a real palindrome node, not a root.
    bool real_node(int node) const {
        return node >= 2 && node < (int)tree.size();
    }

    // Solves: Count distinct palindromic substrings.
    int distinct() const {
        return (int)tree.size() - 2;
    }

    // Solves: Count distinct odd-length palindromic substrings.
    int distinct_odd() const {
        int ans = 0;

        for (int i = 2; i < (int)tree.size(); i++) {
            ans += tree[i].len % 2 == 1;
        }

        return ans;
    }

    // Solves: Count distinct even-length palindromic substrings.
    int distinct_even() const {
        int ans = 0;

        for (int i = 2; i < (int)tree.size(); i++) {
            ans += tree[i].len % 2 == 0;
        }

        return ans;
    }

    // Solves: Return the first stored range [l, r] of a palindrome node.
    pair<int, int> get_range(int node) const {
        int r = tree[node].first_pos;
        int l = r - tree[node].len + 1;

        return {l, r};
    }

    // Solves: Reconstruct the palindrome represented by a node.
    string get_palindrome(int node) const {
        if (!real_node(node)) {
            return "";
        }

        auto [l, r] = get_range(node);
        return s.substr(l, r - l + 1);
    }

    // Solves: Get all real nodes sorted by length.
    vector<int> nodes_by_length(bool decreasing = true) const {
        vector<int> nodes;

        for (int i = 2; i < (int)tree.size(); i++) {
            nodes.push_back(i);
        }

        sort(nodes.begin(), nodes.end(), [&](int a, int b) {
            if (decreasing) {
                return tree[a].len > tree[b].len;
            }

            return tree[a].len < tree[b].len;
        });

        return nodes;
    }

    // Solves: Calculate final occurrence count of every distinct palindrome.
    void count_occurrences() {
        for (int i = 0; i < (int)tree.size(); i++) {
            tree[i].occ = tree[i].raw_occ;
        }

        vector<int> order = nodes_by_length(true);

        for (int v : order) {
            tree[tree[v].link].occ += tree[v].occ;
        }

        occ_ready = true;
    }

    // Solves: Return final occurrence count of a palindrome node.
    long long occurrences_of_node(int node) {
        if (!occ_ready) {
            count_occurrences();
        }

        return tree[node].occ;
    }

    // Solves: Count all palindromic substrings by positions, not only distinct.
    long long total_palindromic_substrings() {
        if (!occ_ready) {
            count_occurrences();
        }

        long long ans = 0;

        for (int i = 2; i < (int)tree.size(); i++) {
            ans += tree[i].occ;
        }

        return ans;
    }

    // Solves: Return the node of the longest palindrome in the whole string.
    int longest_palindrome_node() const {
        return longest_node;
    }

    // Solves: Return the longest palindrome in the whole string.
    string longest_palindrome() const {
        if (longest_node == -1) {
            return "";
        }

        return get_palindrome(longest_node);
    }

    // Solves: Return the longest odd-length palindrome in the whole string.
    string longest_odd_palindrome() const {
        if (longest_odd_node == -1) {
            return "";
        }

        return get_palindrome(longest_odd_node);
    }

    // Solves: Return the longest even-length palindrome in the whole string.
    string longest_even_palindrome() const {
        if (longest_even_node == -1) {
            return "";
        }

        return get_palindrome(longest_even_node);
    }

    // Solves: Return the longest palindromic suffix node of the current full string.
    int longest_suffix_node() const {
        return last;
    }

    // Solves: Return the longest palindromic suffix of the current full string.
    string longest_suffix() const {
        if (s.empty()) {
            return "";
        }

        return get_palindrome(last);
    }

    // Solves: Return length of the longest palindromic suffix of the current full string.
    int max_suffix_len() const {
        if (s.empty()) {
            return 0;
        }

        return tree[last].len;
    }

    // Solves: Return the longest palindromic suffix node after processing position pos.
    int longest_suffix_at(int pos) const {
        return longest_suffix_node_at[pos];
    }

    // Solves: Return the longest palindromic suffix string after processing position pos.
    string longest_suffix_string_at(int pos) const {
        return get_palindrome(longest_suffix_node_at[pos]);
    }

    // Solves: Count distinct palindromes after building prefix s[0..prefix_len-1].
    int distinct_at_prefix(int prefix_len) const {
        if (prefix_len == 0) {
            return 0;
        }

        return distinct_after_prefix[prefix_len - 1];
    }

    // Solves: Return the node created exactly when adding s[pos], or -1 if no new palindrome was created.
    int created_at(int pos) const {
        return created_node_at[pos];
    }

    // Solves: Check whether adding s[pos] created a new distinct palindrome.
    bool created_new_palindrome_at(int pos) const {
        return created_node_at[pos] != -1;
    }

    // Solves: Return the new palindrome created at position pos, or empty string if none was created.
    string new_palindrome_at(int pos) const {
        int node = created_node_at[pos];

        if (node == -1) {
            return "";
        }

        return get_palindrome(node);
    }

    // Solves: Return all palindromic suffix nodes of a given palindrome node.
    vector<int> suffix_chain(int node) const {
        vector<int> chain;

        while (node > 1) {
            chain.push_back(node);
            node = tree[node].link;
        }

        return chain;
    }

    // Solves: Return all palindromic suffixes of the current full string.
    vector<int> current_suffix_chain() const {
        if (s.empty()) {
            return {};
        }

        return suffix_chain(last);
    }

    // Solves: Return all palindromic suffix strings of the current full string.
    vector<string> current_palindromic_suffixes() const {
        vector<string> res;

        for (int node : current_suffix_chain()) {
            res.push_back(get_palindrome(node));
        }

        return res;
    }

    // Solves: Count palindromic suffixes after processing position pos.
    int palindromic_suffix_count_at(int pos) const {
        int node = longest_suffix_node_at[pos];
        return tree[node].suffix_count;
    }

    // Solves: List all palindromic suffix strings after processing position pos.
    vector<string> palindromic_suffixes_at(int pos) const {
        vector<string> res;

        int node = longest_suffix_node_at[pos];

        while (node > 1) {
            res.push_back(get_palindrome(node));
            node = tree[node].link;
        }

        return res;
    }

    // Solves: Return all distinct palindromes.
    vector<string> all_palindromes(bool increasing_length = true) const {
        vector<int> nodes = nodes_by_length(!increasing_length);

        vector<string> res;

        for (int node : nodes) {
            res.push_back(get_palindrome(node));
        }

        return res;
    }

    // Solves: Return all distinct palindromes with final occurrence counts.
    vector<pair<string, long long>> all_palindromes_with_occurrences() {
        if (!occ_ready) {
            count_occurrences();
        }

        vector<pair<string, long long>> res;

        for (int i = 2; i < (int)tree.size(); i++) {
            res.push_back({get_palindrome(i), tree[i].occ});
        }

        return res;
    }

    // Solves: Find the node of a given palindrome string. O(number of distinct palindromes * length).
    int find_palindrome(const string &p) const {
        for (int i = 2; i < (int)tree.size(); i++) {
            if (tree[i].len == (int)p.size() && get_palindrome(i) == p) {
                return i;
            }
        }

        return -1;
    }

    // Solves: Check whether a palindrome exists as a substring.
    bool contains_palindrome(const string &p) const {
        return find_palindrome(p) != -1;
    }

    // Solves: Count how many times a given palindrome string appears.
    long long occurrences_of_palindrome(const string &p) {
        int node = find_palindrome(p);

        if (node == -1) {
            return 0;
        }

        return occurrences_of_node(node);
    }

    // Solves: Return outgoing transitions from a node as {character, target_node}.
    vector<pair<char, int>> transitions(int node) const {
        vector<pair<char, int>> res;

        for (int c = 0; c < SIG; c++) {
            if (tree[node].nxt[c]) {
                res.push_back({char(BASE + c), tree[node].nxt[c]});
            }
        }

        return res;
    }

    // Solves: Return series-link chain, useful in optimized palindrome DP.
    vector<int> series_chain(int node) const {
        vector<int> res;

        while (node > 1) {
            res.push_back(node);
            node = tree[node].series_link;
        }

        return res;
    }

    // Solves: Count distinct palindromes by their length.
    vector<int> distinct_count_by_length() const {
        vector<int> ans(s.size() + 1, 0);

        for (int i = 2; i < (int)tree.size(); i++) {
            ans[tree[i].len]++;
        }

        return ans;
    }

    // Solves: Count total palindrome occurrences by their length.
    vector<long long> occurrence_count_by_length() {
        if (!occ_ready) {
            count_occurrences();
        }

        vector<long long> ans(s.size() + 1, 0);

        for (int i = 2; i < (int)tree.size(); i++) {
            ans[tree[i].len] += tree[i].occ;
        }

        return ans;
    }

    // Solves: Find palindrome maximizing len * occurrence count.
    pair<string, long long> max_len_times_occurrence() {
        if (!occ_ready) {
            count_occurrences();
        }

        long long best = 0;
        int best_node = -1;

        for (int i = 2; i < (int)tree.size(); i++) {
            long long val = 1LL * tree[i].len * tree[i].occ;

            if (val > best) {
                best = val;
                best_node = i;
            }
        }

        if (best_node == -1) {
            return {"", 0};
        }

        return {get_palindrome(best_node), best};
    }

    // Solves: Find the most frequent palindrome.
    pair<string, long long> most_frequent_palindrome() {
        if (!occ_ready) {
            count_occurrences();
        }

        long long best = 0;
        int best_node = -1;

        for (int i = 2; i < (int)tree.size(); i++) {
            if (tree[i].occ > best) {
                best = tree[i].occ;
                best_node = i;
            }
        }

        if (best_node == -1) {
            return {"", 0};
        }

        return {get_palindrome(best_node), best};
    }

    // Solves: Find the longest palindrome that appears at least k times.
    string longest_palindrome_occurring_at_least(long long k) {
        if (!occ_ready) {
            count_occurrences();
        }

        int best_node = -1;

        for (int i = 2; i < (int)tree.size(); i++) {
            if (tree[i].occ >= k) {
                if (best_node == -1 || tree[i].len > tree[best_node].len) {
                    best_node = i;
                }
            }
        }

        if (best_node == -1) {
            return "";
        }

        return get_palindrome(best_node);
    }

    // Solves: Find the longest palindromic prefix of the full string.
    string longest_palindromic_prefix() const {
        int best_node = -1;

        for (int i = 2; i < (int)tree.size(); i++) {
            auto [l, r] = get_range(i);

            if (l == 0) {
                if (best_node == -1 || tree[i].len > tree[best_node].len) {
                    best_node = i;
                }
            }
        }

        if (best_node == -1) {
            return "";
        }

        return get_palindrome(best_node);
    }

    // Solves: Find the longest palindromic suffix of the full string.
    string longest_palindromic_suffix() const {
        return longest_suffix();
    }

    // Solves: Count distinct palindromes that are prefixes of the full string.
    int count_distinct_palindromic_prefixes() const {
        int ans = 0;

        for (int i = 2; i < (int)tree.size(); i++) {
            auto [l, r] = get_range(i);
            ans += l == 0;
        }

        return ans;
    }

    // Solves: Count distinct palindromes that are suffixes of the full string.
    int count_distinct_palindromic_suffixes() const {
        return (int)current_suffix_chain().size();
    }

    // Solves: Check if the whole string is a palindrome.
    bool whole_string_is_palindrome() const {
        if (s.empty()) {
            return true;
        }

        return tree[last].len == (int)s.size();
    }

    // Solves: Minimum number of palindromic pieces to split the string, simple O(total suffix-chain length).
    int min_pal_partition_slow() const {
        int n = (int)s.size();
        const int INF = 1e9;

        vector<int> dp(n + 1, INF);
        dp[0] = 0;

        for (int i = 1; i <= n; i++) {
            int node = longest_suffix_node_at[i - 1];

            while (node > 1) {
                int len = tree[node].len;
                dp[i] = min(dp[i], dp[i - len] + 1);
                node = tree[node].link;
            }
        }

        return dp[n];
    }

    // Solves: Minimum number of palindromic pieces using series links, faster standard Eertree DP.
    int min_pal_partition() const {
        int n = (int)s.size();
        const int INF = 1e9;

        vector<int> dp(n + 1, INF);
        vector<int> best(tree.size(), INF);

        dp[0] = 0;

        for (int i = 1; i <= n; i++) {
            int node = longest_suffix_node_at[i - 1];

            while (node > 1) {
                int sl = tree[node].series_link;
                int pos = i - (tree[sl].len + tree[node].diff);

                best[node] = dp[pos];

                if (tree[node].diff == tree[tree[node].link].diff) {
                    best[node] = min(best[node], best[tree[node].link]);
                }

                dp[i] = min(dp[i], best[node] + 1);

                node = tree[node].series_link;
            }
        }

        return dp[n];
    }

    // Solves: Minimum number of cuts needed so every part is a palindrome.
    int min_pal_cuts() const {
        if (s.empty()) {
            return 0;
        }

        return min_pal_partition() - 1;
    }

    // Solves: Check if the string can be split into at most k palindromic substrings.
    bool can_partition_into_at_most_k_palindromes(int k) const {
        return min_pal_partition() <= k;
    }

    // Solves: Check if the string can be split into exactly k palindromic substrings.
    bool can_partition_into_exactly_k_palindromes(int k) const {
        int n = (int)s.size();
        int mn = min_pal_partition();

        return mn <= k && k <= n;
    }

    // Solves: Count ways to split the string into palindromic substrings, O(total suffix-chain length).
    long long count_pal_partitions(long long mod) const {
        int n = (int)s.size();

        vector<long long> dp(n + 1, 0);
        dp[0] = 1;

        for (int i = 1; i <= n; i++) {
            int node = longest_suffix_node_at[i - 1];

            while (node > 1) {
                int len = tree[node].len;

                dp[i] += dp[i - len];
                dp[i] %= mod;

                node = tree[node].link;
            }
        }

        return dp[n];
    }

    // Solves: Return dp[i] = minimum palindrome partition count for every prefix.
    vector<int> min_pal_partition_for_all_prefixes() const {
        int n = (int)s.size();
        const int INF = 1e9;

        vector<int> dp(n + 1, INF);
        vector<int> best(tree.size(), INF);

        dp[0] = 0;

        for (int i = 1; i <= n; i++) {
            int node = longest_suffix_node_at[i - 1];

            while (node > 1) {
                int sl = tree[node].series_link;
                int pos = i - (tree[sl].len + tree[node].diff);

                best[node] = dp[pos];

                if (tree[node].diff == tree[tree[node].link].diff) {
                    best[node] = min(best[node], best[tree[node].link]);
                }

                dp[i] = min(dp[i], best[node] + 1);

                node = tree[node].series_link;
            }
        }

        return dp;
    }

    // Solves: Return one optimal palindrome partition of the full string. Slow but simple.
    vector<string> restore_one_min_partition_slow() const {
        int n = (int)s.size();
        const int INF = 1e9;

        vector<int> dp(n + 1, INF);
        vector<int> parent(n + 1, -1);

        dp[0] = 0;

        for (int i = 1; i <= n; i++) {
            int node = longest_suffix_node_at[i - 1];

            while (node > 1) {
                int len = tree[node].len;

                if (dp[i - len] + 1 < dp[i]) {
                    dp[i] = dp[i - len] + 1;
                    parent[i] = i - len;
                }

                node = tree[node].link;
            }
        }

        vector<string> parts;

        int cur = n;

        while (cur > 0) {
            int pre = parent[cur];
            parts.push_back(s.substr(pre, cur - pre));
            cur = pre;
        }

        reverse(parts.begin(), parts.end());

        return parts;
    }

    // Solves: Print internal nodes for debugging.
    void debug_print(bool with_occurrences = true) {
        if (with_occurrences && !occ_ready) {
            count_occurrences();
        }

        cout << "Eertree nodes:\n";

        for (int i = 0; i < (int)tree.size(); i++) {
            cout << "node " << i << ": ";

            if (i == 0) {
                cout << "odd root";
            } else if (i == 1) {
                cout << "even root";
            } else {
                cout << '"' << get_palindrome(i) << '"';
            }

            cout << ", len = " << tree[i].len;
            cout << ", link = " << tree[i].link;
            cout << ", diff = " << tree[i].diff;
            cout << ", series_link = " << tree[i].series_link;
            cout << ", suffix_count = " << tree[i].suffix_count;
            cout << ", raw_occ = " << tree[i].raw_occ;

            if (with_occurrences) {
                cout << ", occ = " << tree[i].occ;
            }

            cout << '\n';
        }
    }
};
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string s;
    cin >> s;

    Eertree et;
    et.build(s);

    cout << "Distinct palindromes: " << et.distinct() << '\n';
    cout << "Distinct odd: " << et.distinct_odd() << '\n';
    cout << "Distinct even: " << et.distinct_even() << '\n';

    cout << "Longest palindrome: " << et.longest_palindrome() << '\n';
    cout << "Longest odd palindrome: " << et.longest_odd_palindrome() << '\n';
    cout << "Longest even palindrome: " << et.longest_even_palindrome() << '\n';

    cout << "Longest palindromic prefix: " << et.longest_palindromic_prefix() << '\n';
    cout << "Longest palindromic suffix: " << et.longest_palindromic_suffix() << '\n';

    et.count_occurrences();

    cout << "Total palindromic substrings: " << et.total_palindromic_substrings() << '\n';

    auto [p1, occ] = et.most_frequent_palindrome();
    cout << "Most frequent palindrome: " << p1 << " " << occ << '\n';

    auto [p2, val] = et.max_len_times_occurrence();
    cout << "Best len * occ: " << p2 << " " << val << '\n';

    cout << "Minimum palindrome partition: " << et.min_pal_partition() << '\n';
    cout << "Minimum palindrome cuts: " << et.min_pal_cuts() << '\n';

    cout << "Ways to split into palindromes mod 1e9+7: ";
    cout << et.count_pal_partitions(1000000007LL) << '\n';

    return 0;
}
