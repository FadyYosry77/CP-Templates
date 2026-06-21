const ll inf = 4e18, mod = 998244353; // 998244353
template<int MOD, int RT>
struct mint {
    static const int mod = MOD;
    static constexpr mint rt() { return RT; }
    int v;

    mint() : v(0) {
    }

    mint(ll _v) : v((_v % MOD)) { if(v < 0) v += MOD; }
    explicit operator int() const { return v; }

    mint &operator+=(mint o) {
        if((v += o.v) >= MOD) v -= MOD;
        return *this;
    }

    mint &operator-=(mint o) {
        if((v -= o.v) < 0) v += MOD;
        return *this;
    }

    mint &operator*=(mint o) {
        v = ((ll) v * o.v % MOD);
        return *this;
    }

    friend mint pow(mint a, ll p) {
        assert(p >= 0);
        return p == 0 ? mint(1) : pow(a * a, p / 2) * (p & 1 ? a : mint(1));
    }

    friend mint inv(mint a) {
        assert(a.v != 0);
        return pow(a, MOD - 2);
    }

    friend mint operator+(mint a, mint b) { return a += b; }
    friend mint operator-(mint a, mint b) { return a -= b; }
    friend mint operator*(mint a, mint b) { return a *= b; }
};

using mi = mint<mod, 5>;

vector<mi> fact, invf, pw, ipw;
mi inv2;

void pre(ll n) {
    fact.resize(n + 1), invf.resize(n + 1);
    pw.resize(n + 1), ipw.resize(n + 1);

    fact[0] = mi(1);
    for(ll i = 1;i <= n;i++)fact[i] = fact[i - 1] * mi(i);

    invf[n] = inv(fact[n]);
    for(ll i = n;i > 0;i--)invf[i - 1] = invf[i] * mi(i);

    inv2 = inv(mi(2));

    pw[0] = ipw[0] = mi(1);

    for(ll i = 1;i <= n;i++) {
        pw[i] = pw[i - 1] * mi(2);
        ipw[i] = ipw[i - 1] * inv2;
    }
}

mi ncr(ll n,ll r) {
    if(r < 0 || r > n)return mi(0);
    return fact[n] * invf[r] * invf[n - r];
}
