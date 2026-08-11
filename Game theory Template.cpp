#include <bits/stdc++.h>
using namespace std;

using ll = long long;
__extension__ typedef unsigned __int128 u128;

/*
GameTheory (normal play unless a function says otherwise)

Contains:
- mex, generic W/L and Sprague-Grundy DP, DAG games, winning moves
- Nim, misere Nim, bounded Nim, Moore Nim, Nimble, Staircase Nim
- subtraction, splitting, Kayles, interval and graph W/L/D games
- score games, Euclid, Wythoff, Silver Dollar, Turning Turtles
- Green Hackenbush (tree / parallel edges / general undirected graph)
- divisor games and the complete CSES game set

Conventions:
- A terminal normal-play position is losing and has Grundy value 0.
- A sum of independent impartial games is losing iff their Grundy XOR is 0.
- Callers must give acyclic transitions to recursive DP helpers.
*/

struct GameTheory {
    struct Move {
        int component = -1, to = -1;
    };

    struct ScoreResult {
        ll difference = 0;
        vector<pair<int, int> > moves;
    };

    struct GraphResult {
        // 0 = lose, 1 = win, 2 = draw.
        vector<int> outcome, distance, winning_to;
    };

    // Smallest nonnegative integer absent from a. O(|a|)
    static int mex(const vector<int>&a) {
        vector<char> seen(a.size() + 1);
        for(int x:a) if(0 <= x && x <= (int) a.size()) seen[x] = 1;
        for(int x = 0;;++x) if(!seen[x]) return x;
    }

    // W[x] = whether x has a move to a losing state. States must be topological. O(E)
    static vector<char> indexedWin(const vector<vector<int> >&moves) {
        int n = moves.size();
        vector<char> win(n);
        for(int x = 0;x < n;++x)
            for(int y:moves[x]) {
                assert(0 <= y && y < x);
                if(!win[y]) {
                    win[x] = 1;
                    break;
                }
            }
        return win;
    }

    // Grundy[x] = mex of reachable Grundy values. States must be topological. O(E)
    static vector<int> indexedGrundy(const vector<vector<int> >&moves) {
        int n = moves.size();
        vector<int> g(n), mark(n + 1, -1);
        for(int x = 0;x < n;++x) {
            for(int y:moves[x]) {
                assert(0 <= y && y < x);
                if(g[y] <= n) mark[g[y]] = x;
            }
            while(g[x] <= n && mark[g[x]] == x) ++g[x];
        }
        return g;
    }

    // Grundy values on a DAG. Throws if the graph contains a directed cycle. O(V+E)
    static vector<int> dagGrundy(const vector<vector<int> >&graph) {
        int n = graph.size();
        vector<int> indeg(n), order;
        for(auto&v:graph)
            for(int to:v) {
                assert(0 <= to && to < n);
                ++indeg[to];
            }
        queue<int> q;
        for(int i = 0;i < n;++i) if(!indeg[i]) q.push(i);
        while(!q.empty()) {
            int v = q.front();
            q.pop();
            order.push_back(v);
            for(int to:graph[v]) if(!--indeg[to]) q.push(to);
        }
        if((int) order.size() != n) throw invalid_argument("dagGrundy: cycle");
        vector<int> g(n), seen(n + 1, -1);
        for(int z = n;z--;) {
            int v = order[z];
            for(int to:graph[v]) if(g[to] <= n) seen[g[to]] = v;
            while(g[v] <= n && seen[g[v]] == v) ++g[v];
        }
        return g;
    }

    static int grundyXor(const vector<int>&componentGrundy) {
        int ans = 0;
        for(int x:componentGrundy) ans ^= x;
        return ans;
    }

    // Returns a component move that makes the total XOR zero, or {-1,-1}. O(E)
    static Move grundyWinningMove(const vector<int>&current,
                                  const vector<vector<int> >&reachableValues) {
        assert(current.size() == reachableValues.size());
        int total = grundyXor(current);
        for(int i = 0;i < (int) current.size();++i)
            for(int to:reachableValues[i])
                if((total ^ current[i] ^ to) == 0) return {i, to};
        return {};
    }

    // Ordinary Nim: remove any positive amount from exactly one heap. O(n)
    static bool nimWin(const vector<ll>&piles) {
        ll x = 0;
        for(ll pile:piles) {
            assert(pile >= 0);
            x ^= pile;
        }
        return x != 0;
    }

    // Returns {heap index,new size}; {-1,-1} means the position is losing. O(n)
    static pair<int, ll> nimWinningMove(const vector<ll>&piles) {
        ll x = 0;
        for(ll pile:piles) {
            assert(pile >= 0);
            x ^= pile;
        }
        for(int i = 0;i < (int) piles.size();++i) {
            ll to = piles[i] ^ x;
            if(to < piles[i]) return {i, to};
        }
        return {-1, -1};
    }

    // Misere Nim: taking the last object loses. O(n)
    static bool misereNimWin(const vector<ll>&piles) {
        ll x = 0;
        int nonzero = 0, big = 0;
        for(ll pile:piles) {
            assert(pile >= 0);
            x ^= pile;
            nonzero += pile > 0;
            big += pile > 1;
        }
        if(!nonzero) return true; // Player to move did not take the last object.
        if(!big) return nonzero % 2 == 0;
        return x != 0;
    }

    // Remove 1..k from one heap: g(x)=x%(k+1). This is not for sparse move sets.
    static bool boundedNimWin(const vector<ll>&piles, ll k) {
        assert(k >= 1);
        ll x = 0;
        for(ll pile:piles) {
            assert(pile >= 0);
            x ^= pile % (k + 1);
        }
        return x != 0;
    }

    // Move one coin to a lower-index square. A coin at index i is a Nim heap i.
    static ll nimbleXor(const vector<ll>&coinCount) {
        ll x = 0;
        for(int i = 0;i < (int) coinCount.size();++i) {
            assert(coinCount[i] >= 0);
            if(coinCount[i] & 1) x ^= i;
        }
        return x;
    }

    // Staircase Nim: zero-based odd indices (one-based even stairs) matter. O(n)
    static ll staircaseXor(const vector<ll>&stairs) {
        ll x = 0;
        for(int i = 1;i < (int) stairs.size();i += 2) {
            assert(stairs[i] >= 0);
            x ^= stairs[i];
        }
        return x;
    }

    // Moore Nim: reduce arbitrary positive amounts in 1..k distinct heaps.
    // This is an outcome theorem, not the full SG function and not exact-k Nim.
    static bool mooreNimAtMostKWin(const vector<ll>&piles, int k) {
        assert(k >= 1);
        for(int bit = 0;bit < 63;++bit) {
            int sum = 0;
            for(ll pile:piles) {
                assert(pile >= 0);
                sum += (pile >> bit) & 1;
            }
            if(sum % (k + 1)) return true;
        }
        return false;
    }

    // Turning Turtles: clear H at j, optionally toggle one coin strictly left of j.
    static int turtlesOptionalGrundy(const string&coins) {
        int x = 0;
        for(int j = 0;j < (int) coins.size();++j) {
            assert(coins[j] == 'H' || coins[j] == 'T');
            if(coins[j] == 'H') x ^= j + 1;
        }
        return x;
    }

    // Turning Turtles: clear H at j>0 and must toggle one coin strictly left.
    static int turtlesMandatoryGrundy(const string&coins) {
        int x = 0;
        for(int j = 0;j < (int) coins.size();++j) {
            assert(coins[j] == 'H' || coins[j] == 'T');
            if(coins[j] == 'H') x ^= j;
        }
        return x;
    }

    // One subtraction heap; answer includes state 0. O(maximum * |moves|)
    static vector<char> subtractionWin(int maximum, vector<int> moves) {
        assert(maximum >= 0);
        sort(moves.begin(), moves.end());
        moves.erase(unique(moves.begin(), moves.end()), moves.end());
        for(int d:moves)
            assert(d > 0);
        vector<char> win(maximum + 1);
        for(int x = 1;x <= maximum;++x)
            for(int d:moves) {
                if(d > x) break;
                if(!win[x - d]) {
                    win[x] = 1;
                    break;
                }
            }
        return win;
    }

    // One subtraction heap SG table; answer includes state 0. O(maximum * |moves|)
    static vector<int> subtractionGrundy(int maximum, vector<int> moves) {
        assert(maximum >= 0);
        sort(moves.begin(), moves.end());
        moves.erase(unique(moves.begin(), moves.end()), moves.end());
        for(int d:moves)
            assert(d > 0);
        vector<int> g(maximum + 1), seen(moves.size() + 1, -1);
        for(int x = 1;x <= maximum;++x) {
            for(int d:moves) if(d <= x && g[x - d] < (int) seen.size()) seen[g[x - d]] = x;
            while(g[x] < (int) seen.size() && seen[g[x]] == x) ++g[x];
        }
        return g;
    }

    static int subtractionPilesXor(const vector<int>&piles, const vector<int>&g) {
        int x = 0;
        for(int pile:piles) {
            assert(0 <= pile && pile < (int)g.size());
            x ^= g[pile];
        }
        return x;
    }

    static vector<int> squareMoves(int maximum) {
        vector<int> moves;
        for(ll x = 1;x * x <= maximum;++x) moves.push_back(x * x);
        return moves;
    }

    // g[x] for a heap split into independent (a,b). callback returns legal splits.
    template<class SplitGenerator>
    static vector<int> splittingGrundy(int maximum, SplitGenerator splits) {
        vector<int> g(maximum + 1);
        for(int x = 1;x <= maximum;++x) {
            vector<int> values;
            for(auto [a,b]:splits(x)) {
                assert(0 <= a && a < x && 0 <= b && b < x);
                values.push_back(g[a] ^ g[b]);
            }
            g[x] = mex(values);
        }
        return g;
    }

    // Split into two nonempty unequal heaps. Discovery only: O(maximum^2).
    static vector<int> grundysGameSmall(int maximum) {
        return splittingGrundy(maximum, [](int x) {
            vector<pair<int, int> > out;
            for(int a = 1;a < x - a;++a) if(a != x - a) out.push_back({a, x - a});
            return out;
        });
    }

    // Local blocking game: anchor i also blocks left/right neighbors. O(n^2*types)
    static vector<int> blockedRowGrundy(int maximum,
                                        const vector<pair<int, int> >&moveTypes) {
        vector<int> g(maximum + 1);
        for(int n = 1;n <= maximum;++n) {
            vector<int> values;
            for(int i = 0;i < n;++i)
                for(auto [leftBlock,rightBlock]:moveTypes) {
                    assert(leftBlock >= 0 && rightBlock >= 0);
                    int left = max(0, i - leftBlock);
                    int right = max(0, n - i - 1 - rightBlock);
                    values.push_back(g[left] ^ g[right]);
                }
            g[n] = mex(values);
        }
        return g;
    }

    // Classical Kayles: remove one pin or two adjacent pins. O(n^2)
    static vector<int> kaylesGrundy(int maximum) {
        vector<int> g(maximum + 1);
        for(int n = 1;n <= maximum;++n) {
            vector<int> values;
            for(int i = 0;i < n;++i) values.push_back(g[i] ^ g[n - i - 1]);
            for(int i = 0;i + 1 < n;++i) values.push_back(g[i] ^ g[n - i - 2]);
            g[n] = mex(values);
        }
        return g;
    }

    // Remove one legal point from [l,r), splitting it into two games. O(n^3)
    template<class Legal>
    static vector<vector<int> > intervalRemoveOneGrundy(int n, Legal legal) {
        vector<vector<int> > g(n + 1, vector<int>(n + 1));
        for(int len = 1;len <= n;++len)
            for(int l = 0, r = len;r <= n;++l, ++r) {
                vector<int> values;
                for(int p = l;p < r;++p)
                    if(legal(l, r, p))
                        values.push_back(g[l][p] ^ g[p + 1][r]);
                g[l][r] = mex(values);
            }
        return g;
    }

    // Directed graph outcomes under normal play; unresolved states are draws. O(V+E)
    static GraphResult graphOutcomes(const vector<vector<int> >&graph) {
        int n = graph.size();
        vector<vector<int> > rev(n);
        vector<int> left(n), outcome(n, 2), distance(n), winningTo(n, -1), loseMax(n);
        for(int v = 0;v < n;++v) {
            left[v] = graph[v].size();
            for(int to:graph[v]) {
                assert(0 <= to && to < n);
                rev[to].push_back(v);
            }
        }
        queue<int> q;
        for(int v = 0;v < n;++v) if(!left[v]) outcome[v] = 0, q.push(v);
        while(!q.empty()) {
            int v = q.front();
            q.pop();
            for(int p:rev[v])
                if(outcome[p] == 2) {
                    if(outcome[v] == 0) {
                        outcome[p] = 1;
                        distance[p] = distance[v] + 1;
                        winningTo[p] = v;
                        q.push(p);
                    } else {
                        loseMax[p] = max(loseMax[p], distance[v]);
                        if(!--left[p]) {
                            outcome[p] = 0;
                            distance[p] = loseMax[p] + 1;
                            q.push(p);
                        }
                    }
                }
        }
        return {outcome, distance, winningTo};
    }

    // Removal Game: take left or right, maximize first score - second score. O(n^2)
    static ScoreResult removalGame(const vector<ll>&a) {
        int n = a.size();
        if(!n) return {};
        vector<vector<ll> > dp(n, vector<ll>(n));
        vector<vector<char> > takeLeft(n, vector<char>(n));
        for(int i = 0;i < n;++i) dp[i][i] = a[i], takeLeft[i][i] = 1;
        for(int len = 2;len <= n;++len)
            for(int l = 0, r = len - 1;r < n;++l, ++r) {
                ll left = a[l] - dp[l + 1][r], right = a[r] - dp[l][r - 1];
                if(left >= right) dp[l][r] = left, takeLeft[l][r] = 1;
                else dp[l][r] = right;
            }
        ScoreResult ans;
        ans.difference = dp[0][n - 1];
        for(int l = 0, r = n - 1;l <= r;) {
            if(takeLeft[l][r]) ans.moves.push_back({l, l}), ++l;
            else ans.moves.push_back({r, r}), --r;
        }
        return ans;
    }

    // Euclid's game: replace larger by larger-q*smaller, q>=1. O(log max(a,b))
    static bool euclidWin(ll a, ll b) {
        assert(a >= 0 && b >= 0);
        bool turn = false;
        while(a && b) {
            if(a < b) swap(a, b);
            if(a / b >= 2 || a % b == 0) return !turn;
            a %= b;
            turn = !turn;
        }
        return turn;
    }

    // Wythoff Nim: remove from one heap or equally from both. Exact integer test.
    static bool wythoffLosing(ll a, ll b) {
        assert(a >= 0 && b >= 0);
        if(a > b) swap(a, b);
        unsigned long long k = b - a;
        // a=floor(k*phi): f(a)<=0 and f(a+1)>0 for f(x)=x^2-k*x-k^2.
        u128 kk = (u128) k * k;
        bool left = (u128) a * a <= kk + (u128) a * k;
        bool right = (u128) (a + 1) * (a + 1) > kk + (u128) (a + 1) * k;
        return left && right;
    }

    // Silver Dollar; sorted distinct zero-based coin positions, moves left without jumping.
    static ll silverDollarXor(const vector<ll>&positions) {
        for(int i = 0;i < (int) positions.size();++i) {
            assert(positions[i] >= 0);
            if(i)
                assert(positions[i-1] < positions[i]);
        }
        ll x = 0;
        int i = positions.size() - 1;
        for(;i >= 1;i -= 2) x ^= positions[i] - positions[i - 1] - 1;
        if(i == 0) x ^= positions[0];
        return x;
    }

    // One token may move from a vertex to a child; tree must be rooted at root. O(n)
    static vector<int> treeTokenGrundy(const vector<vector<int> >&tree, int root = 0) {
        int n = tree.size();
        if(!n) return {};
        vector<int> parent(n, -2), order{root};
        parent[root] = -1;
        for(int z = 0;z < (int) order.size();++z) {
            int v = order[z];
            for(int to:tree[v])
                if(to != parent[v]) {
                    if(parent[to] != -2) throw invalid_argument("treeTokenGrundy: not a tree");
                    parent[to] = v;
                    order.push_back(to);
                }
        }
        if((int) order.size() != n) throw invalid_argument("treeTokenGrundy: disconnected");
        vector<int> g(n);
        for(int z = n;z--;) {
            int v = order[z];
            vector<int> values;
            for(int to:tree[v]) if(parent[to] == v) values.push_back(g[to]);
            g[v] = mex(values);
        }
        return g;
    }

    // Green Hackenbush tree: cut one edge; everything disconnected from root falls. O(n)
    static ll greenHackenbushTree(const vector<vector<int> >&tree, int root = 0) {
        int n = tree.size();
        if(!n) return 0;
        vector<int> parent(n, -2), order{root};
        parent[root] = -1;
        for(int z = 0;z < (int) order.size();++z) {
            int v = order[z];
            for(int to:tree[v])
                if(to != parent[v]) {
                    if(parent[to] != -2) throw invalid_argument("greenHackenbushTree: not a tree");
                    parent[to] = v;
                    order.push_back(to);
                }
        }
        if((int) order.size() != n) throw invalid_argument("greenHackenbushTree: disconnected");
        vector<ll> g(n);
        for(int z = n;z--;) {
            int v = order[z];
            for(int to:tree[v]) if(parent[to] == v) g[v] ^= g[to] + 1;
        }
        return g[root];
    }

    // Edge weight = number of parallel green edges, not a chain length. O(n)
    static ll greenHackenbushParallelTree(
        const vector<vector<pair<int, ll> > >&tree, int root = 0) {
        int n = tree.size();
        if(!n) return 0;
        vector<int> parent(n, -2), order{root};
        vector<ll> weight(n);
        parent[root] = -1;
        for(int z = 0;z < (int) order.size();++z) {
            int v = order[z];
            for(auto [to,w]:tree[v])
                if(to != parent[v]) {
                    assert(w >= 1);
                    if(parent[to] != -2) throw invalid_argument("parallel tree: not a tree");
                    parent[to] = v;
                    weight[to] = w;
                    order.push_back(to);
                }
        }
        if((int) order.size() != n) throw invalid_argument("parallel tree: disconnected");
        vector<ll> g(n);
        for(int z = n;z--;) {
            int v = order[z];
            for(auto [to,w]:tree[v])
                if(parent[to] == v) {
                    ll branch = w == 1 ? g[to] + 1 : (w & 1 ? g[to] ^ 1 : g[to]);
                    g[v] ^= branch;
                }
        }
        return g[root];
    }

    // General undirected Green Hackenbush via bridges and parity. O(V+E)
    static ll greenHackenbushGraph(int n, const vector<pair<int, int> >&edges, int root = 0) {
        if(!n) return 0;
        assert(0 <= root && root < n);
        int m = edges.size();
        vector<vector<pair<int, int> > > adj(n);
        for(int id = 0;id < m;++id) {
            auto [a,b] = edges[id];
            assert(0 <= a && a < n && 0 <= b && b < n);
            adj[a].push_back({b, id});
            adj[b].push_back({a, id});
        }
        vector<int> tin(n, -1), low(n), parent(n, -1), parentEdge(n, -1), it(n), order;
        vector<char> bridge(m);
        int timer = 0;
        for(int start = 0;start < n;++start)
            if(tin[start] == -1) {
                vector<int> st{start};
                parent[start] = start;
                while(!st.empty()) {
                    int v = st.back();
                    if(tin[v] == -1) tin[v] = low[v] = timer++;
                    if(it[v] < (int) adj[v].size()) {
                        auto [to,id] = adj[v][it[v]++];
                        if(id == parentEdge[v]) continue;
                        if(tin[to] == -1) parent[to] = v, parentEdge[to] = id, st.push_back(to);
                        else low[v] = min(low[v], tin[to]);
                    } else {
                        st.pop_back();
                        order.push_back(v);
                        if(v != start) {
                            int p = parent[v];
                            low[p] = min(low[p], low[v]);
                            if(low[v] > tin[p]) bridge[parentEdge[v]] = 1;
                        }
                    }
                }
            }
        vector<int> comp(n, -1);
        int cc = 0;
        for(int s = 0;s < n;++s)
            if(comp[s] == -1) {
                vector<int> st{s};
                comp[s] = cc;
                while(!st.empty()) {
                    int v = st.back();
                    st.pop_back();
                    for(auto [to,id]:adj[v])
                        if(!bridge[id] && comp[to] == -1)
                            comp[to] = cc, st.push_back(to);
                }
                ++cc;
            }
        vector<int> parity(cc);
        vector<vector<int> > tree(cc);
        for(int id = 0;id < m;++id) {
            auto [a,b] = edges[id];
            if(bridge[id]) tree[comp[a]].push_back(comp[b]), tree[comp[b]].push_back(comp[a]);
            else parity[comp[a]] ^= 1;
        }
        int cr = comp[root];
        vector<int> par(cc, -2), ord{cr};
        par[cr] = -1;
        for(int z = 0;z < (int) ord.size();++z)
            for(int to:tree[ord[z]]) if(par[to] == -2) par[to] = ord[z], ord.push_back(to);
        vector<ll> g(cc);
        for(int z = ord.size();z--;) {
            int v = ord[z];
            g[v] = parity[v];
            for(int to:tree[v]) if(par[to] == v) g[v] ^= g[to] + 1;
        }
        return g[cr];
    }

    static vector<int> smallestPrimeFactor(int maximum) {
        vector<int> spf(maximum + 1);
        for(int i = 2;i <= maximum;++i)
            if(!spf[i]) {
                spf[i] = i;
                if((ll) i * i <= maximum)
                    for(ll j = (ll) i * i;j <= maximum;j += i) if(!spf[j]) spf[j] = i;
            }
        return spf;
    }

    static vector<pair<int, int> > factorExponents(int n, const vector<int>&spf) {
        assert(n >= 1 && n < (int)spf.size());
        vector<pair<int, int> > out;
        while(n > 1) {
            int p = spf[n], e = 0;
            while(n % p == 0) n /= p, ++e;
            out.push_back({p, e});
        }
        return out;
    }

    // Move n -> n/p for one prime p|n. SG is parity of total prime exponents.
    static int onePrimeFactorGrundy(int n, const vector<int>&spf) {
        int omega = 0;
        for(auto [p,e]:factorExponents(n, spf)) omega += e;
        return omega & 1;
    }

    // Move divides by any positive power of exactly one prime. SG = XOR exponents.
    static int primePowerDivisionGrundy(int n, const vector<int>&spf) {
        int x = 0;
        for(auto [p,e]:factorExponents(n, spf)) x ^= e;
        return x;
    }

    // Move to any proper positive divisor. O(n log n) total transition generation.
    static vector<int> properDivisorGrundy(int maximum) {
        vector<vector<int> > to(maximum + 1);
        for(int d = 1;d <= maximum;++d)
            for(int x = d + d;x <= maximum;x += d) to[x].push_back(d);
        return indexedGrundy(to);
    }

    // Remove one nonzero base-b digit appearing in x. O(maximum * digits)
    static vector<char> subtractDigitWin(int maximum, int base = 10) {
        assert(maximum >= 0 && base >= 2);
        vector<char> win(maximum + 1);
        for(int x = 1;x <= maximum;++x) {
            for(int y = x;y;y /= base) {
                int d = y % base;
                if(d && !win[x - d]) {
                    win[x] = 1;
                    break;
                }
            }
        }
        return win;
    }

    // CSES Stick Game. O(n*|moves|) time, O(n) memory.
    static vector<char> stickGameStates(int n, const vector<int>&moves) {
        return subtractionWin(n, moves);
    }

    static string stickGameString(int n, const vector<int>&moves) {
        auto win = stickGameStates(n, moves);
        string out;
        for(int x = 1;x <= n;++x) out += win[x] ? 'W' : 'L';
        return out;
    }

    static bool nimGameOneWin(const vector<ll>&piles) { return nimWin(piles); }
    static bool nimGameTwoWin(const vector<ll>&piles) { return boundedNimWin(piles, 3); }
    static bool stairGameWin(const vector<ll>&stairs) { return staircaseXor(stairs) != 0; }

    // CSES Grundy's Game, valid only for 1<=n<=1e6.
    // Exact SG computation below 2000 shows the last losing n is 1222;
    // the accepted CSES cutoff treats every n>=2000 as winning.
    static bool grundysGameCsesWin(int n) {
        assert(1 <= n && n <= 1000000);
        if(n >= 2000) return true;
        static const vector<int> g = grundysGameSmall(1999);
        return g[n] != 0;
    }

    // Select any nonempty subset and decrement each selected heap once.
    // Losing iff every heap is even. O(n)
    static bool anotherGameWin(const vector<ll>&piles) {
        for(ll pile:piles) {
            assert(pile >= 0);
            if(pile & 1) return true;
        }
        return false;
    }

    // g(0)=0; from n move to floor(n/2), floor(n/3), or floor(n/6).
    static vector<int> divide236Grundy(int maximum) {
        vector<vector<int> > to(maximum + 1);
        for(int n = 1;n <= maximum;++n) to[n] = {n / 2, n / 3, n / 6};
        return indexedGrundy(to);
    }
};

#ifndef GAME_THEORY_TEMPLATE_NO_MAIN
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Example 1: one subtraction heap / CSES Stick Game.
    cout << GameTheory::stickGameString(10, {1, 3, 4}) << '\n';

    // Example 2: ordinary Nim and one move that leaves XOR zero.
    vector<ll> piles = {5, 7, 2, 5};
    cout << (GameTheory::nimWin(piles) ? "first" : "second") << '\n';
    auto [heap,newSize] = GameTheory::nimWinningMove(piles);
    if(heap != -1)
        cout << "change heap " << heap << " to " << newSize << '\n';

    // Example 3: build a reusable SG table for several subtraction heaps.
    auto g = GameTheory::subtractionGrundy(100, {1, 3, 4});
    vector<int> severalHeaps = {10, 14, 23};
    cout << (GameTheory::subtractionPilesXor(severalHeaps, g)
                 ? "first"
                 : "second") << '\n';

    // Example 4: CSES Stair Game; vectors use zero-based indexing.
    cout << (GameTheory::stairGameWin({0, 2, 1}) ? "first" : "second") << '\n';
}
#endif
