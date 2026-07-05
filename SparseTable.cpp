/*
    Sparse Table Template Guide
    ---------------------------

    Purpose:
    - Static range queries on an array.
    - No updates after build.
    - Best for operations like:
        min, max, gcd, bitwise and, bitwise or

    Indexing:
    - The array is 0-indexed.
    - query(l,r) uses inclusive range [l,r].
    - query_exclusive(l,r) uses range [l,r), where r is not included.
    - Full array query: query(0,n-1).

    Build:
    - sparse<T,F> st(a,operation);
    - Time complexity: O(n log n)
    - Memory complexity: O(n log n)

    Main query:
    - st.query(l,r)
    - Gets merged value in range [l,r].
    - l and r are inclusive.
    - Time complexity: O(1)
    - Works correctly for idempotent operations:
        min, max, gcd, bitwise and, bitwise or

    Log query:
    - st.query_log(l,r)
    - Gets merged value in range [l,r].
    - l and r are inclusive.
    - Time complexity: O(log n)
    - Can be used for associative operations by splitting into non-overlapping blocks.

    Minimum-specific functions:
    - These assume operation is min(x,y):
        query_index(l,r)
        query_pair(l,r)
        has_at_most(l,r,value)
        has_less_than(l,r,value)
        all_greater_than(l,r,value)
        all_at_least(l,r,value)
        first_at_most(l,value)
        first_less_than(l,value)
        last_at_most(r,value)
        last_less_than(r,value)
        max_right_greater_than(l,value)
        min_left_greater_than(r,value)
        length_while_greater_than(l,value)

    Examples:

    Minimum:
        auto operation=[](int x,int y){
            return min(x,y);
        };
        sparse<int,decltype(operation)> st(a,operation);

    Maximum:
        auto operation=[](int x,int y){
            return max(x,y);
        };
        sparse<int,decltype(operation)> st(a,operation);

    GCD:
        auto operation=[](int x,int y){
            return gcd(x,y);
        };
        sparse<int,decltype(operation)> st(a,operation);

    Common mistakes:
    - If input is 1-indexed, convert before query:
        l--, r--;
    - Do not use this template if the array has updates.
    - Do not use minimum-specific helper functions with max/gcd unless you rewrite their logic.
*/
template<typename T, class F>
struct sparse {
    int Log, n;
    vector<vector<T> > table;
    vector<vector<int> > index_table;
    vector<T> values;
    F merge;

    // Builds sparse table from 0-indexed array arr[0..n-1].
    // Range base: all functions use 0-indexed positions.
    // Time complexity: O(n log n)
    explicit sparse(vector<T> arr,F const&merge) : Log(arr.empty() ? 0 : __lg((int) arr.size()) + 1),
                                                   n((int) arr.size()),
                                                   table(Log),
                                                   index_table(Log),
                                                   values(arr),
                                                   merge(merge) {
        if(!n)return;

        table[0] = arr;
        index_table[0].resize(n);

        for(int i = 0;i < n;i++) {
            index_table[0][i] = i;
        }

        for(int l = 1;l < Log;l++) {
            table[l].resize(n + 1 - (1 << l));
            index_table[l].resize(n + 1 - (1 << l));

            for(int i = 0;i < (int) table[l].size();i++) {
                table[l][i] = merge(
                    table[l - 1][i],
                    table[l - 1][i + (1 << (l - 1))]
                );

                int left_index = index_table[l - 1][i];
                int right_index = index_table[l - 1][i + (1 << (l - 1))];

                if(values[left_index] != values[right_index]) {
                    index_table[l][i] = values[left_index] < values[right_index] ? left_index : right_index;
                } else {
                    index_table[l][i] = min(left_index, right_index);
                }
            }
        }
    }

    // Gets merged value in range [l,r].
    // l and r are inclusive.
    // Works in O(1) for idempotent operations like min, max, gcd.
    // Time complexity: O(1)
    T query(int l,int r) {
        if(l > r)return {};
        int len = __lg(r - l + 1);
        return merge(table[len][l], table[len][r - (1 << len) + 1]);
    }

    // Gets merged value in range [l,r).
    // l is inclusive, r is exclusive.
    // Time complexity: O(1)
    T query_exclusive(int l,int r) {
        return query(l, r - 1);
    }

    // Gets merged value in range [l,r] by splitting into powers of two.
    // l and r are inclusive.
    // Useful for non-idempotent operations if merge is associative.
    // Time complexity: O(log n)
    T query_log(int l,int r) {
        if(l > r)return {};

        T answer{};
        bool has_answer = false;

        for(int power = Log - 1;power >= 0;power--) {
            if(l + (1 << power) - 1 <= r) {
                if(!has_answer) {
                    answer = table[power][l];
                    has_answer = true;
                } else {
                    answer = merge(answer, table[power][l]);
                }

                l += 1 << power;
            }
        }

        return answer;
    }

    // Gets index of the first minimum value in range [l,r].
    // l and r are inclusive.
    // This function assumes merge is minimum and T supports comparison.
    // Time complexity: O(1)
    int query_index(int l,int r) {
        if(l > r)return -1;

        int len = __lg(r - l + 1);
        int left_index = index_table[len][l];
        int right_index = index_table[len][r - (1 << len) + 1];

        if(values[left_index] != values[right_index]) {
            return values[left_index] < values[right_index] ? left_index : right_index;
        }

        return min(left_index, right_index);
    }

    // Gets {minimum value, first index of minimum} in range [l,r].
    // l and r are inclusive.
    // This function assumes merge is minimum and T supports comparison.
    // Time complexity: O(1)
    pair<T,int> query_pair(int l,int r) {
        int index = query_index(l, r);
        if(index == -1)return {};
        return {values[index], index};
    }

    // Gets merged value in prefix range [0,r].
    // r is inclusive.
    // Time complexity: O(1)
    T prefix_query(int r) {
        return query(0, r);
    }

    // Gets merged value in suffix range [l,n-1].
    // l is inclusive.
    // Time complexity: O(1)
    T suffix_query(int l) {
        return query(l, n - 1);
    }

    // Checks if minimum value in range [l,r] is <= value.
    // l and r are inclusive.
    // This function assumes merge is minimum.
    // Time complexity: O(1)
    bool has_at_most(int l,int r, T value) {
        return query(l, r) <= value;
    }

    // Checks if minimum value in range [l,r] is < value.
    // l and r are inclusive.
    // This function assumes merge is minimum.
    // Time complexity: O(1)
    bool has_less_than(int l,int r, T value) {
        return query(l, r) < value;
    }

    // Checks if all values in range [l,r] are > value.
    // l and r are inclusive.
    // This function assumes merge is minimum.
    // Time complexity: O(1)
    bool all_greater_than(int l,int r, T value) {
        return query(l, r) > value;
    }

    // Checks if all values in range [l,r] are >= value.
    // l and r are inclusive.
    // This function assumes merge is minimum.
    // Time complexity: O(1)
    bool all_at_least(int l,int r, T value) {
        return query(l, r) >= value;
    }

    // Finds first position pos >= l such that min([l,pos]) <= value.
    // Range checked is [l,pos], inclusive.
    // Returns -1 if no such position exists.
    // This function assumes merge is minimum.
    // Time complexity: O(log n)
    int first_at_most(int l, T value) {
        if(query(l, n - 1) > value)return -1;

        int low = l, high = n - 1;
        while(low < high) {
            int mid = (low + high) / 2;

            if(query(l, mid) <= value)high = mid;
            else low = mid + 1;
        }

        return low;
    }

    // Finds first position pos >= l such that min([l,pos]) < value.
    // Range checked is [l,pos], inclusive.
    // Returns -1 if no such position exists.
    // This function assumes merge is minimum.
    // Time complexity: O(log n)
    int first_less_than(int l, T value) {
        if(query(l, n - 1) >= value)return -1;

        int low = l, high = n - 1;
        while(low < high) {
            int mid = (low + high) / 2;

            if(query(l, mid) < value)high = mid;
            else low = mid + 1;
        }

        return low;
    }

    // Finds last position pos <= r such that min([pos,r]) <= value.
    // Range checked is [pos,r], inclusive.
    // Returns -1 if no such position exists.
    // This function assumes merge is minimum.
    // Time complexity: O(log n)
    int last_at_most(int r, T value) {
        if(query(0, r) > value)return -1;

        int low = 0, high = r;
        while(low < high) {
            int mid = (low + high + 1) / 2;

            if(query(mid, r) <= value)low = mid;
            else high = mid - 1;
        }

        return low;
    }

    // Finds last position pos <= r such that min([pos,r]) < value.
    // Range checked is [pos,r], inclusive.
    // Returns -1 if no such position exists.
    // This function assumes merge is minimum.
    // Time complexity: O(log n)
    int last_less_than(int r, T value) {
        if(query(0, r) >= value)return -1;

        int low = 0, high = r;
        while(low < high) {
            int mid = (low + high + 1) / 2;

            if(query(mid, r) < value)low = mid;
            else high = mid - 1;
        }

        return low;
    }

    // Finds maximum right >= l such that all values in [l,right] are > value.
    // Range checked is [l,right], inclusive.
    // Returns l-1 if values[l] <= value.
    // This function assumes merge is minimum.
    // Time complexity: O(log n)
    int max_right_greater_than(int l, T value) {
        if(query(l, l) <= value)return l - 1;

        int low = l, high = n - 1;
        while(low < high) {
            int mid = (low + high + 1) / 2;

            if(query(l, mid) > value)low = mid;
            else high = mid - 1;
        }

        return low;
    }

    // Finds minimum left <= r such that all values in [left,r] are > value.
    // Range checked is [left,r], inclusive.
    // Returns r+1 if values[r] <= value.
    // This function assumes merge is minimum.
    // Time complexity: O(log n)
    int min_left_greater_than(int r, T value) {
        if(query(r, r) <= value)return r + 1;

        int low = 0, high = r;
        while(low < high) {
            int mid = (low + high) / 2;

            if(query(mid, r) > value)high = mid;
            else low = mid + 1;
        }

        return low;
    }

    // Compares merged value of [l1,r1] with merged value of [l2,r2].
    // All borders are inclusive.
    // Time complexity: O(1)
    bool same_query(int l1,int r1,int l2,int r2) {
        return query(l1, r1) == query(l2, r2);
    }

    // Returns length of longest segment starting at l where all values are > value.
    // Segment is [l,answer], inclusive.
    // This function assumes merge is minimum.
    // Time complexity: O(log n)
    int length_while_greater_than(int l, T value) {
        return max_right_greater_than(l, value) - l + 1;
    }
};

void solve() {
    int n, q;
    cin >> n >> q;

    vector<int> a(n);
    for(int i = 0;i < n;i++)cin >> a[i];

    sparse<int, function<int(int,int)> > st(a, [](int x,int y) {
        return min(x, y);
    });

    while(q--) {
        int l, r;
        cin >> l >> r;

        cout << st.query(l, r) << '\n';
    }
}
