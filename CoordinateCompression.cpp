struct CoordinateCompression {
    vector<ll> cc;

    void init(vector<ll> &a) {
        cc.push_back(-1e18);
        for (auto i: a) cc.push_back(i);
        sort(cc.begin(), cc.end());
        cc.erase(unique(cc.begin(), cc.end()), cc.end());
        cc.push_back(1e18);
    }

    ll idx(ll x) {
        return lower_bound(cc.begin(), cc.end(), x) - cc.begin();
    }

    ll sz() {
        return cc.size();
    }
};
