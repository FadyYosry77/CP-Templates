const int N = 200005;
int n, q, seg;
ll a[N];
vector<ll> comp;

struct Node {
    ll val;
    Node*l,*r;

    Node(ll x) : val(x), l(nullptr), r(nullptr) {
    }

    Node(Node*l2, Node*rr) {
        l = l2;
        r = rr;
        val = 0;
        if(l)val += l->val;
        if(r)val += r->val;
    }

    Node(Node*cp) : val(cp->val), l(cp->l), r(cp->r) {
    }
};

Node*roots[N];

Node*build(int l = 1, int r = seg) {
    if(l == r)return new Node(0LL);

    int mid = (l + r) / 2;

    return new Node(
        build(l, mid),
        build(mid + 1, r)
    );
}

Node*update(Node*node, int val, int pos, int l = 1, int r = seg) {
    if(l == r) {
        return new Node(node->val + val);
    }

    int mid = (l + r) / 2;

    if(pos > mid) {
        return new Node(
            node->l,
            update(node->r, val, pos, mid + 1, r)
        );
    } else {
        return new Node(
            update(node->l, val, pos, l, mid),
            node->r
        );
    }
}

ll query(Node*node, int a, int b, int l = 1, int r = seg) {
    if(!node)return 0;
    if(l > b || r < a)return 0;
    if(l >= a && r <= b)return node->val;

    int mid = (l + r) / 2;

    return query(node->l, a, b, l, mid)
           + query(node->r, a, b, mid + 1, r);
}

int kth(Node*u, Node*v, Node*lc, Node*plc, int k, int l = 1, int r = seg) {
    if(l == r)return l;

    int mid = (l + r) / 2;

    ll leftCount =
            u->l->val +
            v->l->val -
            lc->l->val -
            plc->l->val;

    if(k <= leftCount) {
        return kth(
            u->l,
            v->l,
            lc->l,
            plc->l,
            k,
            l,
            mid
        );
    } else {
        return kth(
            u->r,
            v->r,
            lc->r,
            plc->r,
            k - leftCount,
            mid + 1,
            r
        );
    }
}

// solves: exact compressed index of value x
int get_id(ll x) {
    return lower_bound(all(comp), x) - comp.begin() + 1;
}

// solves: largest compressed index with value <= x
int get_leq_id(ll x) {
    return upper_bound(all(comp), x) - comp.begin();
}

// solves: largest compressed index with value < x
int get_less_id(ll x) {
    return lower_bound(all(comp), x) - comp.begin();
}

// solves: coordinate compression for values a[1..n]
void compress_values() {
    comp.clear();

    for(int i = 1;i <= n;i++) {
        comp.pb(a[i]);
    }

    sort(all(comp));
    comp.erase(unique(all(comp)), comp.end());

    seg = comp.size();
}

// solves: static prefix frequency versions
// use for kth smallest, count <= x, count value frequency in [l,r]
void build_prefix_versions() {
    roots[0] = build();

    for(int i = 1;i <= n;i++) {
        int id = get_id(a[i]);
        roots[i] = update(roots[i - 1], 1, id);
    }
}

// solves: query from difference of two prefix versions
// counts values with compressed index inside [ql,qr] in subarray [l,r]
ll range_query_ids(int l, int r, int ql, int qr) {
    return query(roots[r], ql, qr) - query(roots[l - 1], ql, qr);
}

// solves: frequency of exact value x in subarray [l,r]
ll count_equal(int l, int r, ll x) {
    int id = get_id(x);

    if(id < 1 || id > seg || comp[id - 1] != x)return 0;

    return range_query_ids(l, r, id, id);
}

// solves: count of values <= x in subarray [l,r]
ll count_leq(int l, int r, ll x) {
    int id = get_leq_id(x);

    if(id == 0)return 0;

    return range_query_ids(l, r, 1, id);
}

// solves: count of values < x in subarray [l,r]
ll count_less(int l, int r, ll x) {
    int id = get_less_id(x);

    if(id == 0)return 0;

    return range_query_ids(l, r, 1, id);
}

// solves: count of values >= x in subarray [l,r]
ll count_geq(int l, int r, ll x) {
    int id = lower_bound(all(comp), x) - comp.begin() + 1;

    if(id > seg)return 0;

    return range_query_ids(l, r, id, seg);
}

// solves: count of values > x in subarray [l,r]
ll count_greater(int l, int r, ll x) {
    int id = upper_bound(all(comp), x) - comp.begin() + 1;

    if(id > seg)return 0;

    return range_query_ids(l, r, id, seg);
}

// solves: count values inside value interval [x,y] in subarray [l,r]
ll count_between(int l, int r, ll x, ll y) {
    if(x > y)return 0;

    int left = lower_bound(all(comp), x) - comp.begin() + 1;
    int right = upper_bound(all(comp), y) - comp.begin();

    if(left > right)return 0;
    if(left > seg || right < 1)return 0;

    left = max<ll>(left, 1);
    right = min(right, seg);

    return range_query_ids(l, r, left, right);
}

// solves: kth smallest compressed id in subarray [l,r]
int kth_range_id(Node*u, Node*v, int k, int l = 1, int r = seg) {
    if(l == r)return l;

    int mid = (l + r) / 2;
    ll leftCount = u->l->val - v->l->val;

    if(k <= leftCount) {
        return kth_range_id(u->l, v->l, k, l, mid);
    } else {
        return kth_range_id(u->r, v->r, k - leftCount, mid + 1, r);
    }
}

// solves: kth smallest value in subarray [l,r]
ll kth_smallest(int l, int r, int k) {
    int id = kth_range_id(roots[r], roots[l - 1], k);
    return comp[id - 1];
}

// solves: kth largest value in subarray [l,r]
ll kth_largest(int l, int r, int k) {
    int len = r - l + 1;
    return kth_smallest(l, r, len - k + 1);
}

// solves: lower median of subarray [l,r]
ll range_median(int l, int r) {
    int len = r - l + 1;
    return kth_smallest(l, r, (len + 1) / 2);
}

// solves: minimum value present in subarray [l,r]
ll min_present(int l, int r) {
    return kth_smallest(l, r, 1);
}

// solves: maximum value present in subarray [l,r]
ll max_present(int l, int r) {
    return kth_smallest(l, r, r - l + 1);
}

// solves: number of inversions in static array
ll count_inversions_static() {
    compress_values();
    roots[0] = build();

    ll inv = 0;

    for(int i = 1;i <= n;i++) {
        int id = get_id(a[i]);

        if(id < seg)inv += query(roots[i - 1], id + 1, seg);

        roots[i] = update(roots[i - 1], 1, id);
    }

    return inv;
}

// solves: static distinct values in range using previous occurrence idea
// root[i] stores prev[1],prev[2],...,prev[i] as prev+1
void build_distinct_prev_versions() {
    seg = n + 1;
    roots[0] = build();

    map<ll, int> last;

    for(int i = 1;i <= n;i++) {
        int prev = last[a[i]];
        roots[i] = update(roots[i - 1], 1, prev + 1);
        last[a[i]] = i;
    }
}

// solves: number of distinct values in subarray [l,r]
ll distinct_count(int l, int r) {
    return query(roots[r], 1, l) - query(roots[l - 1], 1, l);
}

// solves: checks if every value in [l,r] is distinct, static version
bool all_distinct_static(int l, int r) {
    return distinct_count(l, r) == r - l + 1;
}

// solves: count duplicate positions in [l,r]
// example: [1,2,1,3,2] in [1,5] has 2 duplicate positions
ll duplicate_positions_count(int l, int r) {
    return (r - l + 1) - distinct_count(l, r);
}

// solves: static distinct count using latest active occurrence idea
// root[i] keeps only latest occurrence of every value active among prefix 1..i
void build_latest_active_versions() {
    seg = n;
    roots[0] = build();

    map<ll, int> last;

    for(int i = 1;i <= n;i++) {
        roots[i] = roots[i - 1];

        if(last[a[i]]) {
            roots[i] = update(roots[i], -1, last[a[i]]);
        }

        roots[i] = update(roots[i], 1, i);
        last[a[i]] = i;
    }
}

// solves: number of distinct values in [l,r] using latest active versions
ll distinct_count_latest(int l, int r) {
    return query(roots[r], l, r);
}

// solves: persistent array / versioned range sum initial build
void build_versioned_array() {
    seg = n;
    roots[0] = build();

    for(int i = 1;i <= n;i++) {
        roots[0] = update(roots[0], a[i], i);
    }
}

// solves: creates new version by adding delta to position pos from old root
Node*new_point_add_version(Node*old_root, int pos, ll delta) {
    return update(old_root, delta, pos);
}

// solves: creates new version by setting position pos to value val from old root
Node*new_point_set_version(Node*old_root, int pos, ll val) {
    ll old = query(old_root, pos, pos);
    return update(old_root, val - old, pos);
}

// solves: range sum in any version root
ll version_range_sum(Node*root, int l, int r) {
    return query(root, l, r);
}

// solves: point value in any version root
ll version_point_value(Node*root, int pos) {
    return query(root, pos, pos);
}

// solves: kth active element inside one frequency version
// example: find kth one in a binary/frequency persistent segment tree
int kth_one(Node*node, ll k, int l = 1, int r = seg) {
    if(l == r)return l;

    int mid = (l + r) / 2;
    ll left = node->l ? node->l->val : 0;

    if(k <= left)return kth_one(node->l, k, l, mid);
    return kth_one(node->r, k - left, mid + 1, r);
}

// solves: first active position in [ql,qr] inside one version
int first_active(Node*node, int ql, int qr, int l = 1, int r = seg) {
    if(!node || node->val == 0 || r < ql || qr < l)return -1;
    if(l == r)return l;

    int mid = (l + r) / 2;

    int left = first_active(node->l, ql, qr, l, mid);
    if(left != -1)return left;

    return first_active(node->r, ql, qr, mid + 1, r);
}

// solves: last active position in [ql,qr] inside one version
int last_active(Node*node, int ql, int qr, int l = 1, int r = seg) {
    if(!node || node->val == 0 || r < ql || qr < l)return -1;
    if(l == r)return l;

    int mid = (l + r) / 2;

    int right = last_active(node->r, ql, qr, mid + 1, r);
    if(right != -1)return right;

    return last_active(node->l, ql, qr, l, mid);
}

// solves: merge two frequency trees
// useful in subtree / DSU-on-tree style problems
Node*merge_trees(Node*x, Node*y, int l = 1, int r = seg) {
    if(!x)return y;
    if(!y)return x;

    if(l == r)return new Node(x->val + y->val);

    int mid = (l + r) / 2;

    return new Node(
        merge_trees(x->l, y->l, l, mid),
        merge_trees(x->r, y->r, mid + 1, r)
    );
}

// solves: path query count on tree after roots[u] are built as root-to-node versions
ll path_query_ids(Node*u, Node*v, Node*lc, Node*plc, int ql, int qr) {
    return query(u, ql, qr)
           + query(v, ql, qr)
           - query(lc, ql, qr)
           - query(plc, ql, qr);
}

// solves: kth smallest value on tree path when four roots are already known
ll kth_on_path_value(Node*u, Node*v, Node*lc, Node*plc, int k) {
    int id = kth(u, v, lc, plc, k);
    return comp[id - 1];
}

// solves: count values <= x on tree path when four roots are already known
ll count_leq_on_path(Node*u, Node*v, Node*lc, Node*plc, ll x) {
    int id = get_leq_id(x);

    if(id == 0)return 0;

    return path_query_ids(u, v, lc, plc, 1, id);
}

// solves: count values < x on tree path when four roots are already known
ll count_less_on_path(Node*u, Node*v, Node*lc, Node*plc, ll x) {
    int id = get_less_id(x);

    if(id == 0)return 0;

    return path_query_ids(u, v, lc, plc, 1, id);
}

// solves: count values >= x on tree path when four roots are already known
ll count_geq_on_path(Node*u, Node*v, Node*lc, Node*plc, ll x) {
    int id = lower_bound(all(comp), x) - comp.begin() + 1;

    if(id > seg)return 0;

    return path_query_ids(u, v, lc, plc, id, seg);
}

// solves: count values > x on tree path when four roots are already known
ll count_greater_on_path(Node*u, Node*v, Node*lc, Node*plc, ll x) {
    int id = upper_bound(all(comp), x) - comp.begin() + 1;

    if(id > seg)return 0;

    return path_query_ids(u, v, lc, plc, id, seg);
}

// solves: count exact value x on tree path when four roots are already known
ll count_equal_on_path(Node*u, Node*v, Node*lc, Node*plc, ll x) {
    int id = get_id(x);

    if(id < 1 || id > seg || comp[id - 1] != x)return 0;

    return path_query_ids(u, v, lc, plc, id, id);
}

void solve() {
    int mode;
    cin >> mode;

    // mode 1:
    // kth smallest in subarray [l,r]
    // input per query: l r k
    // output: kth smallest value in a[l..r]
    if(mode == 1) {
        cin >> n >> q;

        for(int i = 1;i <= n;i++) {
            cin >> a[i];
        }

        compress_values();

        roots[0] = build();

        for(int i = 1;i <= n;i++) {
            int id = get_id(a[i]);
            roots[i] = update(roots[i - 1], 1, id);
        }

        while(q--) {
            int l, r, k;
            cin >> l >> r >> k;

            int id = kth_range(roots[r], roots[l - 1], k);

            cout << comp[id - 1] << '\n';
        }
    }

    // mode 2:
    // count values <= x in subarray [l,r]
    // input per query: l r x
    // output: number of elements <= x in a[l..r]
    else if(mode == 2) {
        cin >> n >> q;

        for(int i = 1;i <= n;i++) {
            cin >> a[i];
        }

        compress_values();

        roots[0] = build();

        for(int i = 1;i <= n;i++) {
            int id = get_id(a[i]);
            roots[i] = update(roots[i - 1], 1, id);
        }

        while(q--) {
            int l, r;
            ll x;
            cin >> l >> r >> x;

            int id = get_leq_id(x);

            if(id == 0)cout << 0 << '\n';
            else cout << query(roots[r], 1, id) - query(roots[l - 1], 1, id) << '\n';
        }
    }

    // mode 3:
    // count values by comparison with x in subarray [l,r]
    // input per query: type l r x
    // type 1: count values < x
    // type 2: count values <= x
    // type 3: count values > x
    // type 4: count values >= x
    else if(mode == 3) {
        cin >> n >> q;

        for(int i = 1;i <= n;i++) {
            cin >> a[i];
        }

        compress_values();

        roots[0] = build();

        for(int i = 1;i <= n;i++) {
            int id = get_id(a[i]);
            roots[i] = update(roots[i - 1], 1, id);
        }

        while(q--) {
            int type, l, r;
            ll x;
            cin >> type >> l >> r >> x;

            if(type == 1) {
                int id = get_less_id(x);
                if(id == 0)cout << 0 << '\n';
                else cout << query(roots[r], 1, id) - query(roots[l - 1], 1, id) << '\n';
            } else if(type == 2) {
                int id = get_leq_id(x);
                if(id == 0)cout << 0 << '\n';
                else cout << query(roots[r], 1, id) - query(roots[l - 1], 1, id) << '\n';
            } else if(type == 3) {
                int id = upper_bound(comp.begin(), comp.end(), x) - comp.begin() + 1;
                if(id > seg)cout << 0 << '\n';
                else cout << query(roots[r], id, seg) - query(roots[l - 1], id, seg) << '\n';
            } else {
                int id = lower_bound(comp.begin(), comp.end(), x) - comp.begin() + 1;
                if(id > seg)cout << 0 << '\n';
                else cout << query(roots[r], id, seg) - query(roots[l - 1], id, seg) << '\n';
            }
        }
    }

    // mode 4:
    // count values inside value interval [x,y] in subarray [l,r]
    // input per query: l r x y
    // output: number of elements x <= a[i] <= y inside [l,r]
    else if(mode == 4) {
        cin >> n >> q;

        for(int i = 1;i <= n;i++) {
            cin >> a[i];
        }

        compress_values();

        roots[0] = build();

        for(int i = 1;i <= n;i++) {
            int id = get_id(a[i]);
            roots[i] = update(roots[i - 1], 1, id);
        }

        while(q--) {
            int l, r;
            ll x, y;
            cin >> l >> r >> x >> y;

            if(x > y) {
                cout << 0 << '\n';
                continue;
            }

            int left = lower_bound(comp.begin(), comp.end(), x) - comp.begin() + 1;
            int right = upper_bound(comp.begin(), comp.end(), y) - comp.begin();

            if(left > right || left > seg || right < 1) {
                cout << 0 << '\n';
                continue;
            }

            left = max(left, 1);
            right = min(right, seg);

            cout << query(roots[r], left, right) - query(roots[l - 1], left, right) << '\n';
        }
    }

    // mode 5:
    // frequency of exact value x in subarray [l,r]
    // input per query: l r x
    // output: count of x inside a[l..r]
    else if(mode == 5) {
        cin >> n >> q;

        for(int i = 1;i <= n;i++) {
            cin >> a[i];
        }

        compress_values();

        roots[0] = build();

        for(int i = 1;i <= n;i++) {
            int id = get_id(a[i]);
            roots[i] = update(roots[i - 1], 1, id);
        }

        while(q--) {
            int l, r;
            ll x;
            cin >> l >> r >> x;

            int id = get_id(x);

            if(id < 1 || id > seg || comp[id - 1] != x) {
                cout << 0 << '\n';
            } else {
                cout << query(roots[r], id, id) - query(roots[l - 1], id, id) << '\n';
            }
        }
    }

    // mode 6:
    // order-statistics helpers in subarray [l,r]
    // input per query: type l r k
    // type 1: kth smallest, uses k
    // type 2: kth largest, uses k
    // type 3: median, ignores k
    // type 4: minimum, ignores k
    // type 5: maximum, ignores k
    else if(mode == 6) {
        cin >> n >> q;

        for(int i = 1;i <= n;i++) {
            cin >> a[i];
        }

        compress_values();

        roots[0] = build();

        for(int i = 1;i <= n;i++) {
            int id = get_id(a[i]);
            roots[i] = update(roots[i - 1], 1, id);
        }

        while(q--) {
            int type, l, r, k;
            cin >> type >> l >> r >> k;

            int len = r - l + 1;
            int need;

            if(type == 1)need = k;
            else if(type == 2)need = len - k + 1;
            else if(type == 3)need = (len + 1) / 2;
            else if(type == 4)need = 1;
            else need = len;

            int id = kth_range(roots[r], roots[l - 1], need);

            cout << comp[id - 1] << '\n';
        }
    }

    // mode 7:
    // distinct values in static subarray [l,r] using previous occurrence
    // input per query: l r
    // output: number of distinct values in a[l..r]
    else if(mode == 7) {
        cin >> n >> q;

        for(int i = 1;i <= n;i++) {
            cin >> a[i];
        }

        seg = n + 1;
        roots[0] = build();

        map<ll, int> last;

        for(int i = 1;i <= n;i++) {
            int prev = last[a[i]];
            roots[i] = update(roots[i - 1], 1, prev + 1);
            last[a[i]] = i;
        }

        while(q--) {
            int l, r;
            cin >> l >> r;

            cout << query(roots[r], 1, l) - query(roots[l - 1], 1, l) << '\n';
        }
    }

    // mode 8:
    // static distinctness check in subarray [l,r]
    // input per query: l r
    // output: YES if all values in a[l..r] are distinct, otherwise NO
    else if(mode == 8) {
        cin >> n >> q;

        for(int i = 1;i <= n;i++) {
            cin >> a[i];
        }

        seg = n + 1;
        roots[0] = build();

        map<ll, int> last;

        for(int i = 1;i <= n;i++) {
            int prev = last[a[i]];
            roots[i] = update(roots[i - 1], 1, prev + 1);
            last[a[i]] = i;
        }

        while(q--) {
            int l, r;
            cin >> l >> r;

            ll distinct = query(roots[r], 1, l) - query(roots[l - 1], 1, l);

            cout << (distinct == r - l + 1 ? "YES" : "NO") << '\n';
        }
    }

    // mode 9:
    // distinct values in static subarray [l,r] using latest active occurrence
    // input per query: l r
    // output: number of distinct values in a[l..r]
    else if(mode == 9) {
        cin >> n >> q;

        for(int i = 1;i <= n;i++) {
            cin >> a[i];
        }

        seg = n;
        roots[0] = build();

        map<ll, int> last;

        for(int i = 1;i <= n;i++) {
            roots[i] = roots[i - 1];

            if(last[a[i]]) {
                roots[i] = update(roots[i], -1, last[a[i]]);
            }

            roots[i] = update(roots[i], 1, i);
            last[a[i]] = i;
        }

        while(q--) {
            int l, r;
            cin >> l >> r;

            cout << query(roots[r], l, r) << '\n';
        }
    }

    // mode 10:
    // persistent versioned array with point assignment and range sum
    // initial version is 0
    // input per query:
    // type 1: 1 ver pos val -> create a new version from version ver after setting a[pos]=val
    // type 2: 2 ver l r     -> query sum in version ver over [l,r]
    else if(mode == 10) {
        cin >> n >> q;

        seg = n;
        roots[0] = build();

        for(int i = 1;i <= n;i++) {
            ll x;
            cin >> x;
            roots[0] = update(roots[0], x, i);
        }

        int versions = 0;

        while(q--) {
            int type;
            cin >> type;

            if(type == 1) {
                int ver, pos;
                ll val;
                cin >> ver >> pos >> val;

                versions++;
                roots[versions] = set_value(roots[ver], val, pos);
            } else {
                int ver, l, r;
                cin >> ver >> l >> r;

                cout << query(roots[ver], l, r) << '\n';
            }
        }
    }

    // mode 11:
    // persistent versioned array with point add and range sum
    // initial version is 0
    // input per query:
    // type 1: 1 ver pos delta -> create a new version from version ver after adding delta to a[pos]
    // type 2: 2 ver l r       -> query sum in version ver over [l,r]
    else if(mode == 11) {
        cin >> n >> q;

        seg = n;
        roots[0] = build();

        for(int i = 1;i <= n;i++) {
            ll x;
            cin >> x;
            roots[0] = update(roots[0], x, i);
        }

        int versions = 0;

        while(q--) {
            int type;
            cin >> type;

            if(type == 1) {
                int ver, pos;
                ll delta;
                cin >> ver >> pos >> delta;

                versions++;
                roots[versions] = update(roots[ver], delta, pos);
            } else {
                int ver, l, r;
                cin >> ver >> l >> r;

                cout << query(roots[ver], l, r) << '\n';
            }
        }
    }

    // mode 12:
    // persistent binary/frequency array: first or last active position in range
    // input per query:
    // type 1: 1 ver pos val -> create new version from ver after setting a[pos]=val
    // type 2: 2 ver l r     -> first active position in [l,r]
    // type 3: 3 ver l r     -> last active position in [l,r]
    else if(mode == 12) {
        cin >> n >> q;

        seg = n;
        roots[0] = build();

        for(int i = 1;i <= n;i++) {
            ll x;
            cin >> x;
            roots[0] = update(roots[0], x, i);
        }

        int versions = 0;

        while(q--) {
            int type;
            cin >> type;

            if(type == 1) {
                int ver, pos;
                ll val;
                cin >> ver >> pos >> val;

                versions++;
                roots[versions] = set_value(roots[ver], val, pos);
            } else if(type == 2) {
                int ver, l, r;
                cin >> ver >> l >> r;

                cout << first_active(roots[ver], l, r) << '\n';
            } else {
                int ver, l, r;
                cin >> ver >> l >> r;

                cout << last_active(roots[ver], l, r) << '\n';
            }
        }
    }

    // mode 13:
    // persistent frequency/binary array: kth active element
    // input per query:
    // type 1: 1 ver pos delta -> create new version from ver after adding delta at pos
    // type 2: 2 ver k         -> find kth active/frequency position in version ver
    else if(mode == 13) {
        cin >> n >> q;

        seg = n;
        roots[0] = build();

        for(int i = 1;i <= n;i++) {
            ll x;
            cin >> x;
            roots[0] = update(roots[0], x, i);
        }

        int versions = 0;

        while(q--) {
            int type;
            cin >> type;

            if(type == 1) {
                int ver, pos;
                ll delta;
                cin >> ver >> pos >> delta;

                versions++;
                roots[versions] = update(roots[ver], delta, pos);
            } else {
                int ver;
                ll k;
                cin >> ver >> k;

                cout << kth_one(roots[ver], k) << '\n';
            }
        }
    }

    // mode 14:
    // count inversions in a static array
    // input:
    // n
    // array
    // output: number of pairs (i,j) such that i<j and a[i]>a[j]
    else if(mode == 14) {
        cin >> n;

        for(int i = 1;i <= n;i++) {
            cin >> a[i];
        }

        compress_values();

        roots[0] = build();

        ll inv = 0;

        for(int i = 1;i <= n;i++) {
            int id = get_id(a[i]);

            if(id < seg)inv += query(roots[i - 1], id + 1, seg);

            roots[i] = update(roots[i - 1], 1, id);
        }

        cout << inv << '\n';
    }

    // mode 15:
    // kth smallest value on a tree path u-v
    // input:
    // n q
    // node values
    // n-1 edges
    // query: u v k
    // output: kth smallest value on path u-v
    else if(mode == 15) {
        cin >> n >> q;

        for(int i = 1;i <= n;i++) {
            cin >> a[i];
        }

        compress_values();

        vector<vector<int> > g(n + 1);

        for(int i = 1;i < n;i++) {
            int u, v;
            cin >> u >> v;

            g[u].push_back(v);
            g[v].push_back(u);
        }

        int lg = 1;
        while((1 << lg) <= n)lg++;

        vector<int> dep(n + 1, 0);
        vector<vector<int> > up(n + 1, vector<int>(lg, 0));

        roots[0] = build();

        function<void(int, int)> dfs = [&](int u, int p) {
            up[u][0] = p;
            dep[u] = dep[p] + 1;

            for(int i = 1;i < lg;i++) {
                up[u][i] = up[up[u][i - 1]][i - 1];
            }

            int id = get_id(a[u]);
            roots[u] = update(roots[p], 1, id);

            for(auto v:g[u]) {
                if(v == p)continue;
                dfs(v, u);
            }
        };

        auto jump = [&](int u, int k) {
            for(int i = 0;i < lg;i++) {
                if(k & (1 << i))u = up[u][i];
            }

            return u;
        };

        auto lca = [&](int u, int v) {
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
        };

        dfs(1, 0);

        while(q--) {
            int u, v, k;
            cin >> u >> v >> k;

            int lc = lca(u, v);
            int plc = up[lc][0];

            int id = kth(roots[u], roots[v], roots[lc], roots[plc], k);

            cout << comp[id - 1] << '\n';
        }
    }
}
