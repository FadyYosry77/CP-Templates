/*
NUMBER THEORY REFERENCE -- standalone C++17/C++20 printable file

CONTAINS
- GCD, checked LCM, extended GCD, Diophantine equations.
- Modular arithmetic, inverses, general CRT, and linear congruences.
- Eratosthenes, linear sieve, SPF, Miller-Rabin, Pollard Rho, segmented sieve.
- Divisor functions, Mobius/phi transforms, coprime and exact-GCD counting.
- Harmonic floor sums, exact roots, square-free and perfect-power methods.
- Subarray/subsequence GCD and LCM state compression.
- Fibonacci, Carmichael lambda, multiplicative order, primitive roots.
- Extended BSGS discrete logarithm and Tonelli-Shanks modular square roots.
- Comprehensive divisor/factorization toolkit and frequency transforms.

NOT INCLUDED
- Combination-specific material is in combinatorics_template.cpp.

FACTORING CHOICE
- Many bounded values: LinearSieve + factorize_spf.
- A few moderate values: factorize or factorize(n,primes).
- Arbitrary signed-64-bit values: factorize64 / factorize_64.

BUILD
    g++ -std=c++17 -O2 number_theory_template.cpp
*/

#include <bits/stdc++.h>
using namespace std;

using ll = long long;
__extension__ typedef __int128 i128;
__extension__ typedef unsigned __int128 u128;


// Problem:
//     Compute gcd of two signed integers.
// Mathematical idea:
//     gcd(a,b)=gcd(b,a%b); every remainder preserves common divisors.
// State / Precomputation:
//     None.
// Returns:
//     Nonnegative gcd, with gcd(0,0)=0.
// Complexity:
//     O(log min(|a|,|b|)).
// Assumptions:
//     Absolute result fits ll; LLONG_MIN is excluded.
// Notes:
//     Every integer linear combination ax+by is divisible by this gcd.
ll nt_gcd(ll a, ll b) {
	a = abs(a); b = abs(b);
	while (b) { ll r = a%b; a = b; b = r; }
	return a;
}

ll nt_gcd(const vector<ll>&a) {
	ll g = 0; for (ll x : a) g = nt_gcd(g,x); return g;
}

// Problem:
//     Compute lcm without silently overflowing.
// Mathematical idea:
//     lcm(a,b)=|a/gcd(a,b)*b|; divide before multiplying.
// State / Precomputation:
//     None.
// Returns:
//     {fits,value}; lcm(0,a)=0.
// Complexity:
//     O(log min(|a|,|b|)).
// Assumptions:
//     Inputs exclude LLONG_MIN.
// Notes:
//     Dividing first reduces but does not eliminate overflow, hence i128 check.
pair<bool,ll> nt_lcm_checked(ll a, ll b) {
	if (!a || !b) return {true,0};
	i128 value = (i128)(a/nt_gcd(a,b))*b; if (value < 0) value = -value;
	if (value > LLONG_MAX) return {false,0};
	return {true,(ll)value};
}

// Problem:
//     Mathematical floor/ceiling of a/b for signed values.
// Mathematical idea:
//     Correct C++ truncation when remainder and denominator have opposite/equal signs.
// State / Precomputation:
//     None.
// Returns:
//     floor(a/b) or ceil(a/b).
// Complexity:
//     O(1).
// Assumptions:
//     b!=0 and result fits ll.
// Notes:
//     Essential for signed Diophantine range bounds.
ll floor_div(ll a, ll b) {
	assert(b != 0); ll q = a/b, r = a%b;
	if (r && ((r > 0) != (b > 0))) q--;
	return q;
}

ll ceil_div(ll a, ll b) {
	assert(b != 0); ll q = a/b, r = a%b;
	if (r && ((r > 0) == (b > 0))) q++;
	return q;
}

i128 floor_div128(i128 a, i128 b) {
	assert(b != 0); i128 q = a/b, r = a%b;
	if (r && ((r > 0) != (b > 0))) q--;
	return q;
}

i128 ceil_div128(i128 a, i128 b) {
	assert(b != 0); i128 q = a/b, r = a%b;
	if (r && ((r > 0) == (b > 0))) q++;
	return q;
}

// Problem:
//     Find Bezout coefficients for ax+by=gcd(a,b).
// Mathematical idea:
//     Reverse Euclid's remainder substitutions.
// State / Precomputation:
//     None.
// Returns:
//     Nonnegative gcd and writes one pair (x,y).
// Complexity:
//     O(log min(|a|,|b|)).
// Assumptions:
//     Coefficients fit ll.
// Notes:
//     For a=b=0 returns g=0,x=1,y=0.
ll extgcd(ll a, ll b, ll&x, ll&y) {
	ll aa = abs(a), bb = abs(b); i128 oldR = aa, r = bb, oldX = 1, xx = 0, oldY = 0, yy = 1;
	while (r) {
		i128 q = oldR/r;
		tie(oldR,r) = pair<i128,i128>{r,oldR-q*r};
		tie(oldX,xx) = pair<i128,i128>{xx,oldX-q*xx};
		tie(oldY,yy) = pair<i128,i128>{yy,oldY-q*yy};
	}
	x = (ll)(a < 0 ? -oldX : oldX); y = (ll)(b < 0 ? -oldY : oldY);
	return (ll)oldR;
}

// Problem:
//     Find one integer solution of ax+by=c.
// Mathematical idea:
//     Solutions exist iff gcd(a,b) divides c; scale a Bezout identity.
// State / Precomputation:
//     None.
// Returns:
//     Whether solvable; writes x,y,g.
// Complexity:
//     O(log min(|a|,|b|)).
// Assumptions:
//     A returned solution fits ll.
// Notes:
//     All solutions are x=x0+k*b/g, y=y0-k*a/g.
bool find_any_solution(ll a, ll b, ll c, ll&x, ll&y, ll&g) {
	if (!a && !b) { g = 0; x = y = 0; return c == 0; }
	g = extgcd(a,b,x,y); if (c%g) return false;
	i128 scale = c/g, xx = (i128)x*scale, yy = (i128)y*scale;
	if (xx < LLONG_MIN || xx > LLONG_MAX || yy < LLONG_MIN || yy > LLONG_MAX) return false;
	x = xx; y = yy; return true;
}

// Problem:
//     Count solutions of ax+by=c inside inclusive x/y ranges.
// Mathematical idea:
//     Intersect the integer-k intervals induced by x=x0+k*b/g and y=y0-k*a/g.
// State / Precomputation:
//     One Bezout solution.
// Returns:
//     Number of solutions, saturated at LLONG_MAX.
// Complexity:
//     O(log min(|a|,|b|)).
// Assumptions:
//     lx<=rx and ly<=ry.
// Notes:
//     Handles zero and negative coefficients explicitly.
ll count_diophantine_solutions(ll a, ll b, ll c,
		ll lx, ll rx, ll ly, ll ry) {
	if (lx > rx || ly > ry) return 0;
	if (!a && !b) {
		if (c) return 0;
		i128 count = ((i128)rx-lx+1)*((i128)ry-ly+1);
		return count > LLONG_MAX ? LLONG_MAX : (ll)count;
	}
	if (!a) {
		if (c%b || c/b < ly || c/b > ry) return 0;
		i128 count = (i128)rx-lx+1; return count > LLONG_MAX ? LLONG_MAX : (ll)count;
	}
	if (!b) {
		if (c%a || c/a < lx || c/a > rx) return 0;
		i128 count = (i128)ry-ly+1; return count > LLONG_MAX ? LLONG_MAX : (ll)count;
	}
	ll x0,y0,g; if (!find_any_solution(a,b,c,x0,y0,g)) return 0;
	i128 low = -((i128)1<<126), high = ((i128)1<<126);
	auto restrict = [&](i128 base, i128 step, i128 left, i128 right) {
		if (step > 0) low = max(low,ceil_div128(left-base,step)), high = min(high,floor_div128(right-base,step));
		else low = max(low,ceil_div128(right-base,step)), high = min(high,floor_div128(left-base,step));
	};
	restrict(x0,b/g,lx,rx); restrict(y0,-a/g,ly,ry);
	if (low > high) return 0;
	i128 count = high-low+1;
	return count > LLONG_MAX ? LLONG_MAX : (ll)count;
}


ll mod_normalize(ll x, ll mod) {
	assert(mod > 0); x %= mod; if (x < 0) x += mod; return x;
}

// Problem:
//     Multiply modulo mod without signed 64-bit overflow.
// Mathematical idea:
//     Evaluate the product in 128-bit arithmetic.
// State / Precomputation:
//     None.
// Returns:
//     a*b mod mod in [0,mod).
// Complexity:
//     O(1) on common compilers.
// Assumptions:
//     0<mod<=LLONG_MAX.
// Notes:
//     Required by 64-bit Miller-Rabin and Pollard Rho.
ll mul_mod(ll a, ll b, ll mod) {
	return (i128)mod_normalize(a,mod)*mod_normalize(b,mod)%mod;
}

ll power_mod(ll a, unsigned long long exponent, ll mod) {
	assert(mod > 0); ll result = 1%mod; a = mod_normalize(a,mod);
	while (exponent) {
		if (exponent&1) result = mul_mod(result,a,mod);
		a = mul_mod(a,a,mod); exponent >>= 1;
	}
	return result;
}

// Problem:
//     Invert a modulo an arbitrary positive modulus.
// Mathematical idea:
//     ax+mod*y=1 is solvable exactly when gcd(a,mod)=1.
// State / Precomputation:
//     Extended gcd.
// Returns:
//     Inverse in [0,mod), or -1 when it does not exist.
// Complexity:
//     O(log mod).
// Assumptions:
//     mod>1.
// Notes:
//     Unlike a^(mod-2), this does not assume mod is prime.
ll inverse_mod(ll a, ll mod) {
	assert(mod > 1); ll x,y,g = extgcd(a,mod,x,y);
	return g == 1 ? mod_normalize(x,mod) : -1;
}

vector<ll> inverses_prime(int n, ll primeMod) {
	assert(n >= 0 && n < primeMod); vector<ll> inverse(n+1);
	if (n) inverse[1] = 1;
	for (int i = 2; i <= n; i++) inverse[i] = primeMod-mul_mod(primeMod/i,inverse[primeMod%i],primeMod);
	return inverse;
}

// Problem:
//     Compute a^n and 1+a+...+a^(n-1) modulo mod.
// Mathematical idea:
//     Double blocks; no division by a-1 is needed.
// State / Precomputation:
//     Recursive exponent-halving identities.
// Returns:
//     {a^n,geometric sum of n terms}.
// Complexity:
//     O(log n).
// Assumptions:
//     mod>0.
// Notes:
//     Safe even when a-1 has no modular inverse.
pair<ll,ll> power_geometric_sum(ll a, unsigned long long n, ll mod) {
	if (!n) return {1%mod,0};
	auto [power,sum] = power_geometric_sum(a,n/2,mod);
	ll doubledPower = mul_mod(power,power,mod);
	ll doubledSum = mul_mod(sum,(1+power)%mod,mod);
	if (!(n&1)) return {doubledPower,doubledSum};
	return {mul_mod(doubledPower,a,mod),(doubledSum+doubledPower)%mod};
}

// Problem:
//     Merge x=a (mod m) and x=b (mod n), allowing non-coprime moduli.
// Mathematical idea:
//     Consistency requires a=b modulo gcd(m,n); solve the reduced linear congruence.
// State / Precomputation:
//     Extended gcd.
// Returns:
//     {smallest nonnegative solution,lcm}; second=-1 on inconsistency/ll overflow.
// Complexity:
//     O(log min(m,n)).
// Assumptions:
//     m,n>0.
// Notes:
//     The solution is unique modulo lcm(m,n).
pair<ll,ll> crt_merge(ll a, ll m, ll b, ll n) {
	assert(m > 0 && n > 0); a = mod_normalize(a,m); b = mod_normalize(b,n);
	ll x,y,g = extgcd(m,n,x,y); ll difference = b-a;
	if (difference%g) return {0,-1};
	i128 modulus = (i128)(m/g)*n; if (modulus > LLONG_MAX) return {0,-1};
	ll reduced = n/g;
	ll k = mul_mod(difference/g,mod_normalize(x,reduced),reduced);
	ll result = (ll)(((i128)a+(i128)m*k)%modulus);
	return {result,(ll)modulus};
}

pair<ll,ll> crt_system(const vector<ll>&remainders, const vector<ll>&moduli) {
	assert(remainders.size() == moduli.size()); pair<ll,ll> answer{0,1};
	for (int i = 0; i < (int)remainders.size() && answer.second != -1; i++)
		answer = crt_merge(answer.first,answer.second,remainders[i],moduli[i]);
	return answer;
}

// Problem:
//     Solve a*x=b (mod modulus).
// Mathematical idea:
//     Divide by g=gcd(a,m); the reduced coefficient is invertible.
// State / Precomputation:
//     Extended-gcd inverse.
// Returns:
//     {one solution, reduced modulus}; second=-1 if impossible.
// Complexity:
//     O(log modulus).
// Assumptions:
//     modulus>0.
// Notes:
//     If solvable there are g incongruent solutions modulo the original modulus.
pair<ll,ll> solve_linear_congruence(ll a, ll b, ll modulus) {
	assert(modulus > 0); ll g = nt_gcd(a,modulus); if (b%g) return {0,-1};
	ll reduced = modulus/g; if (reduced == 1) return {0,1};
	ll inverse = inverse_mod(mod_normalize(a/g,reduced),reduced);
	return {mul_mod(b/g,inverse,reduced),reduced};
}


// Problem:
//     List all primes and primality flags up to n.
// Mathematical idea:
//     Mark multiples starting at p^2.
// State / Precomputation:
//     Boolean primality array.
// Returns:
//     {primes,isPrime-as-char}.
// Complexity:
//     O(n log log n) time, O(n) memory.
// Assumptions:
//     n>=0.
// Notes:
//     Prefer this simple sieve when SPF/phi/mu are unnecessary.
pair<vector<int>,vector<char>> eratosthenes(int n) {
	vector<char> isPrime(n+1,true); if (n >= 0) isPrime[0] = false; if (n >= 1) isPrime[1] = false;
	for (ll p = 2; p*p <= n; p++) if (isPrime[p])
		for (ll x = p*p; x <= n; x += p) isPrime[x] = false;
	vector<int> primes; for (int x = 2; x <= n; x++) if (isPrime[x]) primes.push_back(x);
	return {primes,isPrime};
}

struct LinearSieve {
	int n;
	vector<int> primes,spf,phi,mu,omega,bigOmega;

	// Problem:
	//     Precompute primes, SPF, phi, mu and prime-factor counts for every x<=n.
	// Mathematical idea:
	//     Generate each composite once using its smallest prime factor.
	// State / Precomputation:
	//     spf[x], phi[x], mu[x], omega distinct and bigOmega with multiplicity.
	// Returns:
	//     Data members of this struct.
	// Complexity:
	//     O(n) time and O(n) memory.
	// Assumptions:
	//     n>=1.
	// Notes:
	//     phi(1)=mu(1)=1.
	LinearSieve(int limit = 1) : n(max(1,limit)), spf(n+1), phi(n+1),
		mu(n+1), omega(n+1), bigOmega(n+1) {
		spf[1] = phi[1] = mu[1] = 1;
		for (int x = 2; x <= n; x++) {
			if (!spf[x]) {
				spf[x] = x; phi[x] = x-1; mu[x] = -1; omega[x] = bigOmega[x] = 1; primes.push_back(x);
			}
			for (int p : primes) {
				ll value = (ll)x*p; if (value > n) break; spf[value] = p; bigOmega[value] = bigOmega[x]+1;
				if (x%p == 0) { phi[value] = phi[x]*p; mu[value] = 0; omega[value] = omega[x]; break; }
				phi[value] = phi[x]*(p-1); mu[value] = -mu[x]; omega[value] = omega[x]+1;
			}
		}
	}
};

// Purpose:
//     Prime-factorize a bounded integer using its smallest-prime-factor table.
// Input:
//     Positive x and an SPF array covering x.
// Returns:
//     Sorted {prime,exponent} pairs.
// Complexity:
//     O(number of factors with multiplicity), at most O(log x).
// Assumptions:
//     1<=x<spf.size().
// Notes:
//     Preferred when many inputs share a known upper bound.
vector<pair<int,int>> factorize_spf(int x, const vector<int>&spf) {
	assert(1 <= x && x < (int)spf.size()); vector<pair<int,int>> factors;
	while (x > 1) {
		int p = spf[x], exponent = 0;
		while (x%p == 0) x /= p, exponent++;
		factors.push_back({p,exponent});
	}
	return factors;
}

// Problem:
//     Factor one moderate integer using primes up to its square root.
// Mathematical idea:
//     Every composite remainder has a prime factor <=sqrt(remainder).
// State / Precomputation:
//     Optional prime list; pass sieve primes through sqrt(n).
// Returns:
//     Sorted {prime,exponent} pairs.
// Complexity:
//     O(pi(sqrt n)) trial divisions in the worst case.
// Assumptions:
//     n>=1 and primes cover sqrt(n).
// Notes:
//     Usually enough for a few values near 1e12.
vector<pair<ll,int>> factorize_trial(ll n, const vector<int>&primes) {
	assert(n >= 1); vector<pair<ll,int>> factors;
	for (ll p : primes) {
		if (p > n/p) break;
		if (n%p) continue;
		int exponent = 0;
		while (n%p == 0) n /= p, exponent++;
		factors.push_back({p,exponent});
	}
	if (n > 1) factors.push_back({n,1});
	return factors;
}

// Problem:
//     Test primality of any positive signed-64-bit integer.
// Mathematical idea:
//     Strong probable-prime tests with a proven deterministic 64-bit base set.
// State / Precomputation:
//     Bases {2,325,9375,28178,450775,9780504,1795265022}.
// Returns:
//     True exactly for primes.
// Complexity:
//     O(7*log n) modular squarings.
// Assumptions:
//     n fits signed ll.
// Notes:
//     mul_mod uses i128; 1 and negative values are composite.
bool is_prime_64(ll n) {
	if (n < 2) return false;
	for (ll p : {2LL,3LL,5LL,7LL,11LL,13LL,17LL,19LL,23LL,29LL,31LL,37LL}) {
		if (n%p == 0) return n == p;
	}
	unsigned long long d = n-1; int s = 0;
	while (!(d&1)) d >>= 1, s++;
	for (ll base : {2LL,325LL,9375LL,28178LL,450775LL,9780504LL,1795265022LL}) {
		if (base%n == 0) continue;
		ll x = power_mod(base,d,n);
		if (x == 1 || x == n-1) continue;
		bool composite = true;
		for (int r = 1; r < s; r++) {
			x = mul_mod(x,x,n);
			if (x == n-1) { composite = false; break; }
		}
		if (composite) return false;
	}
	return true;
}

mt19937_64 pollard_rng(chrono::steady_clock::now().time_since_epoch().count());

ll pollard_rho_factor(ll n) {
	if (n%2 == 0) return 2;
	if (n%3 == 0) return 3;
	while (true) {
		ll c = uniform_int_distribution<ll>(1,n-1)(pollard_rng);
		ll x = uniform_int_distribution<ll>(0,n-1)(pollard_rng), y = x, divisor = 1;
		auto next = [&](ll value) { return (ll)(((i128)mul_mod(value,value,n)+c)%n); };
		while (divisor == 1) {
			x = next(x); y = next(next(y));
			divisor = nt_gcd(x >= y ? x-y : y-x,n);
		}
		if (divisor != n) return divisor;
	}
}

void factorize_64_raw(ll n, vector<ll>&factors) {
	if (n == 1) return;
	if (is_prime_64(n)) { factors.push_back(n); return; }
	ll divisor = pollard_rho_factor(n);
	factorize_64_raw(divisor,factors); factorize_64_raw(n/divisor,factors);
}

// Problem:
//     Factor an arbitrary positive signed-64-bit integer.
// Mathematical idea:
//     Miller-Rabin recognizes prime leaves; Pollard Rho finds nontrivial factors.
// State / Precomputation:
//     Randomized polynomial iteration.
// Returns:
//     Sorted {prime,exponent} pairs.
// Complexity:
//     Randomized expected sublinear factor search; no simple strict contest bound.
// Assumptions:
//     1<=n<=LLONG_MAX.
// Notes:
//     Trial division is simpler and often faster when n<=1e12 or there are few queries.
vector<pair<ll,int>> factorize_64(ll n) {
	assert(n >= 1); vector<ll> raw; factorize_64_raw(n,raw); sort(raw.begin(),raw.end());
	vector<pair<ll,int>> factors;
	for (ll p : raw) {
		if (factors.empty() || factors.back().first != p) factors.push_back({p,0});
		factors.back().second++;
	}
	return factors;
}

// Problem:
//     Enumerate primes in a manageable segment [left,right] with large coordinates.
// Mathematical idea:
//     Mark the segment using base primes through sqrt(right).
// State / Precomputation:
//     A temporary boolean segment and an ordinary base sieve.
// Returns:
//     All primes in the interval.
// Complexity:#include "number_theory_template.cpp"
//     Roughly O(sqrt(R) log log R + (R-L+1) log log R).
// Assumptions:
//     0<=left<=right and the segment fits memory.
// Notes:
//     Start at max(p^2,ceil(left/p)*p), so p itself is not marked composite.
vector<ll> segmented_primes(ll left, ll right) {
	assert(0 <= left && left <= right); ll root = sqrtl((long double)right);
	while ((i128)(root+1)*(root+1) <= right) root++;
	while ((i128)root*root > right) root--;
	assert(root <= INT_MAX);
	auto [primes,flags] = eratosthenes((int)root); (void)flags;
	vector<char> isPrime(right-left+1,true);
	for (ll p : primes) {
		i128 first = max((i128)p*p,ceil_div128(left,p)*(i128)p);
		for (i128 x = first; x <= right; x += p) isPrime[(ll)x-left] = false;
	}
	for (ll x = left; x <= min(right,1LL); x++) isPrime[x-left] = false;
	vector<ll> answer; for (ll x = left; x <= right; x++) if (isPrime[x-left]) answer.push_back(x);
	return answer;
}


ll divisor_count(const vector<pair<ll,int>>&factors) {
	i128 answer = 1;
	for (auto [prime,exponent] : factors) { (void)prime; answer *= exponent+1; }
	assert(answer <= LLONG_MAX); return answer;
}

// Problem:
//     Compute exact sum of positive divisors from a prime factorization.
// Mathematical idea:
//     sigma(product p^e)=product(1+p+...+p^e).
// State / Precomputation:
//     Prime exponents.
// Returns:
//     {fits,value}.
// Complexity:
//     O(sum exponents).
// Assumptions:
//     Factorization represents a positive integer.
// Notes:
//     i128 detects overflow.
pair<bool,ll> divisor_sum_checked(const vector<pair<ll,int>>&factors) {
	i128 answer = 1;
	for (auto [prime,exponent] : factors) {
		i128 term = 1, power = 1;
		for (int i = 0; i < exponent; i++) power *= prime, term += power;
		answer *= term; if (answer > LLONG_MAX) return {false,0};
	}
	return {true,(ll)answer};
}

ll divisor_sum_mod(const vector<pair<ll,int>>&factors, ll mod) {
	ll answer = 1%mod;
	for (auto [prime,exponent] : factors)
		answer = mul_mod(answer,power_geometric_sum(prime,exponent+1,mod).second,mod);
	return answer;
}

// Problem:
//     Generate all positive divisors of a factored integer.
// Mathematical idea:
//     Choose exponent 0..e independently for each prime.
// State / Precomputation:
//     Growing product list.
// Returns:
//     Sorted divisor vector.
// Complexity:
//     O(tau(n) log tau(n)) including sorting.
// Assumptions:
//     Every generated divisor fits ll.
// Notes:
//     Divisor-lattice DP should compress these divisors instead of allocating up to huge n.
vector<ll> divisors_from_factors(const vector<pair<ll,int>>&factors) {
	vector<ll> divisors{1};
	for (auto [prime,exponent] : factors) {
		int old = divisors.size(); ll power = 1;
		for (int e = 1; e <= exponent; e++) {
			power *= prime;
			for (int i = 0; i < old; i++) divisors.push_back(divisors[i]*power);
		}
	}
	sort(divisors.begin(),divisors.end()); return divisors;
}

ll euler_phi_from_factors(ll n, const vector<pair<ll,int>>&factors) {
	for (auto [prime,exponent] : factors) { (void)exponent; n -= n/prime; }
	return n;
}

ll euler_phi(ll n) {
	assert(n >= 1); return euler_phi_from_factors(n,factorize_64(n));
}

vector<ll> divisor_count_sieve(int n) {
	vector<ll> count(n+1);
	for (int divisor = 1; divisor <= n; divisor++)
		for (int multiple = divisor; multiple <= n; multiple += divisor) count[multiple]++;
	return count;
}

vector<ll> divisor_sum_sieve(int n) {
	vector<ll> sum(n+1);
	for (int divisor = 1; divisor <= n; divisor++)
		for (int multiple = divisor; multiple <= n; multiple += divisor) sum[multiple] += divisor;
	return sum;
}

vector<ll> divisor_zeta_transform(vector<ll> values) {
	int n = values.size()-1; vector<ll> original = values; fill(values.begin(),values.end(),0);
	for (int divisor = 1; divisor <= n; divisor++)
		for (int multiple = divisor; multiple <= n; multiple += divisor)
			values[multiple] += original[divisor];
	return values;
}

vector<ll> divisor_mobius_transform(vector<ll> sums) {
	int n = sums.size()-1;
	for (int divisor = 1; divisor <= n; divisor++)
		for (int multiple = divisor+divisor; multiple <= n; multiple += divisor)
			sums[multiple] -= sums[divisor];
	return sums;
}

vector<ll> count_input_multiples(const vector<int>&values) {
	int maximum = values.empty() ? 0 : *max_element(values.begin(),values.end()); vector<ll> frequency(maximum+1);
	for (int x : values) { assert(x >= 1); frequency[x]++; }
	vector<ll> count(maximum+1);
	for (int divisor = 1; divisor <= maximum; divisor++)
		for (int multiple = divisor; multiple <= maximum; multiple += divisor) count[divisor] += frequency[multiple];
	return count;
}


// Problem:
//     Count unordered index pairs whose values are coprime.
// Mathematical idea:
//     1[gcd=1]=sum_{d|gcd}mu[d].
// State / Precomputation:
//     Counts of input values divisible by each d and mu[d].
// Returns:
//     Number of i<j with gcd(a[i],a[j])=1.
// Complexity:
//     O(M log M+n), M=max(a).
// Assumptions:
//     All values are positive and bounded enough for arrays.
// Notes:
//     C(count,2) means unordenumber_theory_templatered distinct indices.
ll count_coprime_pairs(const vector<int>&values) {
	if (values.empty()) return 0;
	int maximum = *max_element(values.begin(),values.end());
	LinearSieve sieve(maximum); auto multiples = count_input_multiples(values); ll answer = 0;
	for (int d = 1; d <= maximum; d++) answer += (ll)sieve.mu[d]*multiples[d]*(multiples[d]-1)/2;
	return answer;
}

// Problem:
//     Count unordered pairs for every exact gcd.
// Mathematical idea:
//     Pairs divisible by d include exact gcd d,2d,3d,...; subtract multiples descending.
// State / Precomputation:
//     Number of input values divisible by every d.
// Returns:
//     exact[g]=number of i<j with gcd(a[i],a[j])=g.
// Complexity:
//     O(M log M+n).
// Assumptions:
//     Positive bounded values.
// Notes:
//     This subtract-multiples inversion is often clearer than explicit Mobius.
vector<ll> pair_gcd_frequency(const vector<int>&values) {
	if (values.empty()) return {0};
	int maximum = *max_element(values.begin(),values.end());
	auto multiples = count_input_multiples(values); vector<ll> exact(maximum+1);
	for (int d = maximum; d >= 1; d--) {
		exact[d] = multiples[d]*(multiples[d]-1)/2;
		for (int multiple = d+d; multiple <= maximum; multiple += d) exact[d] -= exact[multiple];
	}
	return exact;
}

// Problem:
//     Sum gcd(a[i],a[j]) over all unordered pairs i<j.
// Mathematical idea:
//     gcd(x,y)=sum_{d|gcd(x,y)}phi(d).
// State / Precomputation:
//     Multiples counts and phi values.
// Returns:
//     Pair gcd sum.
// Complexity:
//     O(M log M+n).
// Assumptions:
//     Positive bounded values and answer fits ll.
// Notes:
//     For ordered i!=j multiply by 2; including i=j additionally adds sum a[i].
ll sum_pair_gcd(const vector<int>&values) {
	if (values.empty()) return 0;
	int maximum = *max_element(values.begin(),values.end());
	LinearSieve sieve(maximum); auto multiples = count_input_multiples(values); i128 answer = 0;
	for (int d = 1; d <= maximum; d++) answer += (i128)sieve.phi[d]*multiples[d]*(multiples[d]-1)/2;
	assert(answer <= LLONG_MAX); return answer;
}

// Problem:
//     Count ordered pairs (a,b), 1<=a,b<=n, with gcd(a,b)=1.
// Mathematical idea:
//     Pairs with a<b correspond to phi(b); reflect across the diagonal and include (1,1).
// State / Precomputation:
//     phi[1..n].
// Returns:
//     Number of ordered coprime lattice points in the n by n square.
// Complexity:
//     O(n) after a sieve, O(n) sieve construction.
// Assumptions:
//     n>=1.
// Notes:
//     Answer is 1+2*sum_{x=2}^n phi(x).
ll coprime_lattice_points(int n) {
	assert(n >= 1); LinearSieve sieve(n); i128 answer = 1;
	for (int x = 2; x <= n; x++) answer += 2LL*sieve.phi[x];
	assert(answer <= LLONG_MAX); return answer;
}

// Problem:
//     Count ordered positive pairs with gcd=G and lcm=L.
// Mathematical idea:
//     a=G*x,b=G*y, gcd(x,y)=1 and x*y=L/G; each prime power chooses one side.
// State / Precomputation:
//     Factorization of L/G.
// Returns:
//     2^omega(L/G), or 0 when G does not divide L.
// Complexity:
//     Factorization cost.
// Assumptions:
//     G,L>0.
// Notes:
//     Divide by 2 for unordered distinct-side pairs when the problem requests that convention.
ll count_pairs_with_gcd_lcm(ll gcdValue, ll lcmValue) {
	assert(gcdValue > 0 && lcmValue > 0); if (lcmValue%gcdValue) return 0;
	int distinct = factorize_64(lcmValue/gcdValue).size(); return 1LL<<distinct;
}

// Problem:
//     Count integers in [1,limit] coprime to fixed n.
// Mathematical idea:
//     Inclusion-exclusion over the distinct prime factors of n.
// State / Precomputation:
//     Distinct primes of n.
// Returns:
//     Count of x with gcd(x,n)=1.
// Complexity:
//     O(2^omega(n)) after factorization.
// Assumptions:
//     n>=1,limit>=0.
// Notes:
//     For one fixed n this is simpler than a global Mobius sieve.
ll count_coprime_to_fixed(ll limit, ll n) {
	assert(limit >= 0 && n >= 1); vector<ll> primes;
	for (auto [prime,exponent] : factorize_64(n)) { (void)exponent; primes.push_back(prime); }
	ll notCoprime = 0; int count = primes.size();
	for (int mask = 1; mask < (1<<count); mask++) {
		i128 product = 1; int bits = 0;
		for (int i = 0; i < count; i++) if (mask>>i&1) product *= primes[i], bits++;
		ll term = limit/(ll)product; notCoprime += bits&1 ? term : -term;
	}
	return limit-notCoprime;
}

// Problem:
//     Count integers <=limit divisible by at least one supplied positive number.
// Mathematical idea:
//     Subset inclusion-exclusion; intersection means divisibility by subset lcm.
// State / Precomputation:
//     None beyond subset enumeration.
// Returns:
//     Union size.
// Complexity:
//     O(2^k*k).
// Assumptions:
//     k is small (normally <=20), limit>=0.
// Notes:
//     LCM values above limit are stopped before overflow.
ll count_divisible_by_any(ll limit, vector<ll>numbers) {
	assert(limit >= 0); numbers.erase(remove_if(numbers.begin(),numbers.end(),
		[](ll x) { return x <= 0; }),numbers.end());
	sort(numbers.begin(),numbers.end()); numbers.erase(unique(numbers.begin(),numbers.end()),numbers.end());
	int k = numbers.size(); assert(k <= 25); i128 answer = 0;
	for (unsigned long long mask = 1; mask < (1ULL<<k); mask++) {
		ll lcm = 1; int bits = 0; bool tooLarge = false;
		for (int i = 0; i < k; i++) if (mask>>i&1) {
			bits++; ll g = nt_gcd(lcm,numbers[i]); i128 next = (i128)(lcm/g)*numbers[i];
			if (next > limit) { tooLarge = true; break; } lcm = next;
		}
		if (!tooLarge) answer += bits&1 ? limit/lcm : -limit/lcm;
	}
	return answer;
}


// Problem:
//     Compute sum_{i=1}^n floor(n/i).
// Mathematical idea:
//     q=floor(n/l) stays constant through r=floor(n/q).
// State / Precomputation:
//     O(sqrt n) quotient blocks.
// Returns:
//     Divisor summatory function sum_{x<=n}tau(x).
// Complexity:
//     O(sqrt n).
// Assumptions:
//     n>=0 and result fits ll.
// Notes:
//     Both sides count pairs (a,b) with a*b<=n.
ll sum_floor_divisions(ll n) {
	assert(n >= 0); i128 answer = 0;
	for (ll left = 1,right; left <= n; left = right+1) {
		ll quotient = n/left; right = n/quotient; answer += (i128)quotient*(right-left+1);
	}
	assert(answer <= LLONG_MAX); return answer;
}

// Problem:
//     Compute sum_{0<=i<n} floor((a*i+b)/m).
// Mathematical idea:
//     Repeatedly extract whole multiples, then swap axes in the remaining lattice triangle.
// State / Precomputation:
//     Euclidean-style reductions.
// Returns:
//     The floor sum.
// Complexity:
//     O(log max(a,m)).
// Assumptions:
//     n>=0,m>0; a,b may be signed; answer fits ll.
// Notes:
//     Parameter order matches AtCoder: (n,m,a,b).
ll floor_sum(ll n, ll m, ll a, ll b) {
	assert(n >= 0 && m > 0); i128 answer = 0;
	if (a < 0) {
		ll normalized = mod_normalize(a,m);
		answer -= (i128)n*(n-1)/2*(((i128)normalized-a)/m); a = normalized;
	}
	if (b < 0) {
		ll normalized = mod_normalize(b,m);
		answer -= (i128)n*(((i128)normalized-b)/m); b = normalized;
	}
	while (true) {
		if (a >= m) answer += (i128)(n-1)*n*(a/m)/2, a %= m;
		if (b >= m) answer += (i128)n*(b/m), b %= m;
		i128 maximum = (i128)a*n+b; if (maximum < m) break;
		n = maximum/m; b = maximum%m; swap(m,a);
	}
	assert(answer >= LLONG_MIN && answer <= LLONG_MAX); return answer;
}

// Problem:
//     Exact integer square root.
// Mathematical idea:
//     Correct a sqrtl approximation using i128 comparisons.
// State / Precomputation:
//     None.
// Returns:
//     Largest r with r^2<=n.
// Complexity:
//     O(1) corrections in practice.
// Assumptions:
//     n>=0.
// Notes:
//     Never trust floating point without correction.
ll integer_sqrt(ll n) {
	assert(n >= 0); ll root = sqrtl((long double)n);
	while ((i128)(root+1)*(root+1) <= n) root++;
	while ((i128)root*root > n) root--;
	return root;
}

bool power_leq(ll base, int exponent, ll limit) {
	i128 value = 1;
	for (int i = 0; i < exponent; i++) { value *= base; if (value > limit) return false; }
	return true;
}

ll integer_kth_root(ll n, int exponent) {
	assert(n >= 0 && exponent >= 1); if (exponent == 1) return n;
	ll low = 0, high = min(n,3037000500LL)+1;
	while (low+1 < high) {
		ll middle = low+(high-low)/2;
		if (power_leq(middle,exponent,n)) low = middle; else high = middle;
	}
	return low;
}

// Problem:
//     Count square-free positive integers <=n.
// Mathematical idea:
//     1[squarefree(x)]=sum_{d^2|x}mu(d).
// State / Precomputation:
//     mu through floor(sqrt(n)).
// Returns:
//     Square-free count.
// Complexity:
//     O(sqrt n) after O(sqrt n) sieve.
// Assumptions:
//     n>=0 and sqrt(n) fits int memory.
// Notes:
//     Enables binary search for the k-th square-free number.
ll count_squarefree(ll n) {
	if (n <= 0) return 0;
	ll root64 = integer_sqrt(n); assert(root64 <= INT_MAX); int root = root64;
	LinearSieve sieve(root); ll answer = 0;
	for (ll d = 1; d <= root; d++) answer += (ll)sieve.mu[d]*(n/(d*d));
	return answer;
}

ll kth_squarefree(ll k) {
	assert(k >= 1); ll low = 0, high = max(2LL,2*k);
	while (count_squarefree(high) < k) high *= 2;
	while (low+1 < high) {
		ll middle = low+(high-low)/2;
		if (count_squarefree(middle) >= k) high = middle; else low = middle;
	}
	return high;
}

ll squarefree_kernel(ll n) {
	assert(n >= 1); i128 kernel = 1;
	for (auto [prime,exponent] : factorize_64(n)) if (exponent&1) kernel *= prime;
	return kernel;
}

// Problem:
//     Count pairs of positive values whose product is a perfect square.
// Mathematical idea:
//     Products are squares exactly when square-free kernels are equal.
// State / Precomputation:
//     Frequency map of exponent-parity signatures.
// Returns:
//     Number of unordered index pairs.
// Complexity:
//     Factorization cost plus expected O(n) hashing.
// Assumptions:
//     Values are positive; zero/sign handling depends on the problem and is excluded.
// Notes:
//     For bounded values, SPF makes kernel computation much faster.
ll count_square_product_pairs(const vector<ll>&values) {
	unordered_map<ll,ll> frequency; ll answer = 0;
	for (ll value : values) { ll kernel = squarefree_kernel(value); answer += frequency[kernel]++; }
	return answer;
}

// Problem:
//     Find all positive moduli d for which all values are congruent modulo d.
// Mathematical idea:
//     d must divide every difference, hence d divides gcd(|a[i]-a[0]|).
// State / Precomputation:
//     Gcd of differences and its factorization.
// Returns:
//     Sorted divisors; empty means all values equal and every positive modulus works.
// Complexity:
//     O(n log A)+factorization/divisor generation.
// Assumptions:
//     Differences fit ll.
// Notes:
//     Apply additional lower-bound restrictions after receiving the divisor list.
vector<ll> congruence_moduli(const vector<ll>&values) {
	if (values.size() <= 1) return {};
	ll differenceGcd = 0;
	for (ll value : values) differenceGcd = nt_gcd(differenceGcd,value-values[0]);
	if (!differenceGcd) return {};
	return divisors_from_factors(factorize_64(differenceGcd));
}


// Problem:
//     Count occurrences of every gcd among nonempty contiguous subarrays.
// Mathematical idea:
//     Extend gcd states ending at previous index; equal results are merged.
// State / Precomputation:
//     O(log maxA) distinct gcd states per endpoint.
// Returns:
//     Map gcd -> number of subarrays.
// Complexity:
//     Approximately O(n log maxA).
// Assumptions:
//     Values are nonnegative.
// Notes:
//     Gcd strictly falls through divisors, keeping the compressed state list small.
map<ll,ll> subarray_gcd_frequency(const vector<ll>&values) {
	vector<pair<ll,ll>> previous; map<ll,ll> total;
	for (ll value : values) {
		assert(value >= 0); vector<pair<ll,ll>> current{{value,1}};
		for (auto [g,count] : previous) {
			ll next = nt_gcd(g,value);
			if (current.back().first == next) current.back().second += count;
			else current.push_back({next,count});
		}
		for (auto [g,count] : current) total[g] += count;
		previous = move(current);
	}
	return total;
}

ll count_subarrays_gcd(const vector<ll>&values, ll target) {
	auto frequency = subarray_gcd_frequency(values); return frequency[target];
}

// Problem:
//     Count nonempty subsequences by exact gcd for bounded positive values.
// Mathematical idea:
//     Subsequences made only of multiples of d total 2^cnt[d]-1; subtract exact multiples.
// State / Precomputation:
//     Multiples counts and descending divisor inversion.
// Returns:
//     exact[g] modulo mod.
// Complexity:
//     O(M log M+n).
// Assumptions:
//     Positive values, mod>0.
// Notes:
//     This counts index subsequences, including duplicate values distinctly.
vector<ll> count_subsequences_gcd_all(const vector<int>&values, ll mod) {
	assert(mod > 0); if (values.empty()) return {0};
	int maximum = *max_element(values.begin(),values.end()); auto multiples = count_input_multiples(values);
	vector<ll> exact(maximum+1);
	for (int d = maximum; d >= 1; d--) {
		exact[d] = mod_normalize(power_mod(2,multiples[d],mod)-1,mod);
		for (int multiple = d+d; multiple <= maximum; multiple += d)
			exact[d] = mod_normalize(exact[d]-exact[multiple],mod);
	}
	return exact;
}

// Problem:
//     Count subarrays whose lcm equals target.
// Mathematical idea:
//     Values not dividing target reset; merge compressed lcm states capped by target.
// State / Precomputation:
//     Distinct lcm states ending at each position.
// Returns:
//     Number of qualifying subarrays.
// Complexity:
//     O(n*tau(target)*log target) worst-pattern bound, usually much smaller.
// Assumptions:
//     Values and target are positive.
// Notes:
//     Every surviving lcm divides target, bounding the state universe by its divisors.
ll count_subarrays_lcm(const vector<ll>&values, ll target) {
	assert(target > 0); vector<pair<ll,ll>> previous; ll answer = 0;
	for (ll value : values) {
		assert(value > 0);
		if (target%value) { previous.clear(); continue; }
		vector<pair<ll,ll>> current{{value,1}};
		for (auto [old,count] : previous) {
			ll next = old/nt_gcd(old,value)*value;
			if (current.back().first == next) current.back().second += count;
			else current.push_back({next,count});
		}
		for (auto [lcm,count] : current) if (lcm == target) answer += count;
		previous = move(current);
	}
	return answer;
}


ll vp_factorial(ll n, ll prime) {
	assert(n >= 0 && prime >= 2); ll exponent = 0;
	while (n) n /= prime, exponent += n;
	return exponent;
}

ll vp_binomial(ll n, ll r, ll prime) {
	if (r < 0 || r > n) return 0;
	return vp_factorial(n,prime)-vp_factorial(r,prime)-vp_factorial(n-r,prime);
}

ll trailing_zeros_factorial_base(ll n, ll base) {
	assert(n >= 0 && base >= 2); ll answer = LLONG_MAX;
	for (auto [prime,exponent] : factorize_64(base)) answer = min(answer,vp_factorial(n,prime)/exponent);
	return answer;
}

pair<ll,ll> fibonacci(ll n, ll mod) {
	assert(n >= 0 && mod > 0); if (!n) return {0,1%mod};
	auto [a,b] = fibonacci(n/2,mod);
	ll c = mul_mod(a,mod_normalize(2*b-a,mod),mod);
	ll d = (mul_mod(a,a,mod)+mul_mod(b,b,mod))%mod;
	return n&1 ? pair<ll,ll>{d,(c+d)%mod} : pair<ll,ll>{c,d};
}

ll carmichael_lambda(ll n) {
	assert(n >= 1); ll answer = 1;
	for (auto [prime,exponent] : factorize_64(n)) {
		ll value = 1;
		if (prime == 2 && exponent >= 3) value = 1LL<<(exponent-2);
		else {
			for (int i = 1; i < exponent; i++) value *= prime;
			value *= prime-1;
		}
		auto [fits,lcm] = nt_lcm_checked(answer,value); assert(fits); answer = lcm;
	}
	return answer;
}

// Problem:
//     Compute multiplicative order of a modulo m.
// Mathematical idea:
//     The order divides phi(m); remove each prime factor while the power remains 1.
// State / Precomputation:
//     Factorizations of m and phi(m).
// Returns:
//     Smallest positive k with a^k=1 mod m, or -1 if gcd(a,m)!=1.
// Complexity:
//     Factorization plus O(omega(phi)*log phi) modular work.
// Assumptions:
//     m>=2.
// Notes:
//     Carmichael lambda may be a smaller initial multiple, but phi is sufficient.
ll multiplicative_order(ll a, ll modulus) {
	assert(modulus >= 2); if (nt_gcd(a,modulus) != 1) return -1;
	ll order = euler_phi(modulus);
	for (auto [prime,exponent] : factorize_64(order)) {
		(void)exponent;
		while (order%prime == 0 && power_mod(a,order/prime,modulus) == 1) order /= prime;
	}
	return order;
}

// Problem:
//     Find a primitive root modulo a prime.
// Mathematical idea:
//     g has order p-1 iff g^((p-1)/q)!=1 for every prime q|(p-1).
// State / Precomputation:
//     Distinct prime factors of p-1.
// Returns:
//     Smallest primitive root; 1 for p=2.
// Complexity:
//     Factorization plus candidates*omega(p-1)*log p.
// Assumptions:
//     p is prime.
// Notes:
//     Primitive roots for general composite moduli require an existence test.
ll primitive_root_prime(ll primeMod) {
	assert(is_prime_64(primeMod)); if (primeMod == 2) return 1;
	vector<ll> factors; for (auto [prime,exponent] : factorize_64(primeMod-1)) {
		(void)exponent; factors.push_back(prime);
	}
	for (ll root = 2;; root++) {
		bool good = true;
		for (ll prime : factors) if (power_mod(root,(primeMod-1)/prime,primeMod) == 1) { good = false; break; }
		if (good) return root;
	}
}

ll discrete_log_coprime(ll a, ll b, ll modulus) {
	a = mod_normalize(a,modulus); b = mod_normalize(b,modulus);
	if (modulus == 1) return 0;
	ll size = integer_sqrt(modulus)+1;
	unordered_map<ll,ll> baby; ll value = 1;
	for (ll j = 0; j < size; j++) { if (!baby.count(value)) baby[value] = j; value = mul_mod(value,a,modulus); }
	ll giantInverse = inverse_mod(power_mod(a,size,modulus),modulus);
	ll current = b;
	for (ll i = 0; i <= size; i++) {
		auto it = baby.find(current); if (it != baby.end()) return i*size+it->second;
		current = mul_mod(current,giantInverse,modulus);
	}
	return -1;
}

// Problem:
//     Find minimum x>=0 satisfying a^x=b mod m, even when gcd(a,m)!=1.
// Mathematical idea:
//     Strip common gcd factors, then apply baby-step giant-step in the coprime modulus.
// State / Precomputation:
//     Hash table of O(sqrt m) baby powers.
// Returns:
//     Minimum exponent, or -1 if none exists.
// Complexity:
//     O(sqrt m) expected time and memory.
// Assumptions:
//     m>=1.
// Notes:
//     Memory makes BSGS unsuitable for extremely large m near 1e18.
ll discrete_log(ll a, ll b, ll modulus) {
	assert(modulus >= 1); if (modulus == 1) return 0;
	a = mod_normalize(a,modulus); b = mod_normalize(b,modulus);
	ll offset = 0, accumulated = 1%modulus;
	while (true) {
		ll g = nt_gcd(a,modulus); if (g == 1) break;
		if (b == accumulated) return offset;
		if (b%g) return -1;
		b /= g; modulus /= g; offset++;
		accumulated = mul_mod(accumulated,a/g,modulus);
	}
	if (modulus == 1) return offset;
	ll inverse = inverse_mod(accumulated,modulus);
	ll result = discrete_log_coprime(a,mul_mod(b,inverse,modulus),modulus);
	return result == -1 ? -1 : result+offset;
}

// Problem:
//     Solve x^2=value modulo an odd prime.
// Mathematical idea:
//     Tonelli-Shanks reduces the 2-power part of p-1 using a quadratic nonresidue.
// State / Precomputation:
//     Factor p-1=q*2^s.
// Returns:
//     One root in [0,p), or -1 if no root exists.
// Complexity:
//     O(log^2 p) worst-case modular operations, typically near O(log p).
// Assumptions:
//     primeMod is prime.
// Notes:
//     The other nonzero root is p-root; handles p=2 and value=0.
ll mod_sqrt_prime(ll value, ll primeMod) {
	assert(is_prime_64(primeMod)); value = mod_normalize(value,primeMod);
	if (!value || primeMod == 2) return value;
	if (power_mod(value,(primeMod-1)/2,primeMod) != 1) return -1;
	if (primeMod%4 == 3) return power_mod(value,(primeMod+1)/4,primeMod);
	ll q = primeMod-1; int s = 0; while (!(q&1)) q >>= 1, s++;
	ll nonresidue = 2; while (power_mod(nonresidue,(primeMod-1)/2,primeMod) != primeMod-1) nonresidue++;
	ll x = power_mod(value,(q+1)/2,primeMod), t = power_mod(value,q,primeMod);
	ll c = power_mod(nonresidue,q,primeMod); int m = s;
	while (t != 1) {
		int i = 1; ll squared = mul_mod(t,t,primeMod);
		while (i < m && squared != 1) squared = mul_mod(squared,squared,primeMod), i++;
		ll factor = power_mod(c,1ULL<<(m-i-1),primeMod);
		x = mul_mod(x,factor,primeMod); c = mul_mod(factor,factor,primeMod);
		t = mul_mod(t,c,primeMod); m = i;
	}
	return x;
}


// ============================================================================
// DIVISOR / PRIME-FACTORIZATION UTILITY TOOLKIT
//
// small bounded values   -> LinearSieve + factorize_spf
// one/few medium values  -> factorize or factorize(n,primes)
// arbitrary 64-bit value -> factorize64 (Miller-Rabin + Pollard Rho)
// ============================================================================

// Purpose:
//     Prime-factorize one positive integer by direct trial division.
// Input:
//     n>=1.
// Returns:
//     Sorted {prime,exponent} pairs; n=1 gives an empty vector.
// Complexity:
//     O(sqrt n) worst case.
// Assumptions:
//     Intended for one/few moderate values.
// Notes:
//     The condition prime<=n/prime avoids square overflow.
vector<pair<ll,int>> factorize(ll n) {
	assert(n >= 1); vector<pair<ll,int>> factors;
	for (ll prime = 2; prime <= n/prime; prime += prime == 2 ? 1 : 2) {
		if (n%prime) continue;
		int exponent = 0;
		while (n%prime == 0) n /= prime, exponent++;
		factors.push_back({prime,exponent});
	}
	if (n > 1) factors.push_back({n,1});
	return factors;
}

// Purpose:
//     Prime-factorize using an already generated prime list.
// Input:
//     n>=1 and primes covering every prime through sqrt(original n).
// Returns:
//     Sorted {prime,exponent} pairs.
// Complexity:
//     O(pi(sqrt n)) divisibility checks in the worst case.
// Assumptions:
//     The supplied list is sorted, prime, and sufficiently long.
// Notes:
//     Preferred for a few values near 1e12 after one sieve to 1e6.
vector<pair<ll,int>> factorize(ll n, const vector<int>&primes) {
	return factorize_trial(n,primes);
}

// Purpose:
//     Expose arbitrary signed-64-bit factorization under the compact toolkit name.
// Input:
//     n>=1.
// Returns:
//     Sorted {prime,exponent} pairs.
// Complexity:
//     Randomized Pollard-Rho factorization; expected sublinear factor search.
// Assumptions:
//     n fits positive ll.
// Notes:
//     Use only when trial division/SPF is not sufficient.
vector<pair<ll,int>> factorize64(ll n) { return factorize_64(n); }

// Purpose:
//     Expose deterministic signed-64-bit primality with the compact toolkit name.
// Input:
//     Any signed ll value.
// Returns:
//     True exactly when n is prime.
// Complexity:
//     O(log n) modular squarings for a fixed deterministic base set.
// Assumptions:
//     None beyond signed-64-bit range.
// Notes:
//     This aliases the library's Miller-Rabin implementation.
bool is_prime64(ll n) { return is_prime_64(n); }

// Purpose:
//     Factor an arbitrary positive ll into repeated prime factors.
// Input:
//     n>=1.
// Returns:
//     Sorted primes with multiplicity.
// Complexity:
//     Expected Pollard-Rho factorization time plus output size.
// Assumptions:
//     n fits positive ll.
// Notes:
//     Prefer prime_factors_spf for many bounded values.
vector<ll> prime_factors64(ll n) {
	vector<ll> answer;
	for (auto [prime,exponent] : factorize64(n))
		for (int i = 0; i < exponent; i++) answer.push_back(prime);
	return answer;
}

// Purpose:
//     Alias the safe divide-free geometric power/sum primitive.
// Input:
//     Base a, term count n, modulus>0.
// Returns:
//     {a^n,1+a+...+a^(n-1)} modulo modulus.
// Complexity:
//     O(log n).
// Assumptions:
//     Modular products fit the existing mul_mod contract.
// Notes:
//     Works for composite moduli and when a-1 is not invertible.
pair<ll,ll> power_and_sum(ll a, unsigned long long n, ll modulus) {
	return power_geometric_sum(a,n,modulus);
}

// Purpose:
//     Return prime factors with multiplicity.
// Input:
//     Positive n.
// Returns:
//     For example 12 -> {2,2,3}.
// Complexity:
//     O(sqrt n) through factorize.
// Assumptions:
//     n>=1.
// Notes:
//     Use factorize64 first for arbitrary 64-bit inputs.
vector<ll> prime_factors(ll n) {
	vector<ll> answer;
	for (auto [prime,exponent] : factorize(n))
		for (int i = 0; i < exponent; i++) answer.push_back(prime);
	return answer;
}

// Purpose: Factor a bounded integer into repeated primes using SPF.
// Input: n>=1 and an SPF table covering n.
// Returns: Sorted prime factors with multiplicity.
// Complexity: O(bigOmega(n)).
// Assumptions: n<spf.size().
// Notes: This is the repeated-prime companion to factorize_spf.
vector<int> prime_factors_spf(int n, const vector<int>&spf) {
	assert(1 <= n && n < (int)spf.size()); vector<int> answer;
	while (n > 1) answer.push_back(spf[n]), n /= spf[n];
	return answer;
}

// Purpose:
//     Return each distinct prime divisor once.
// Input:
//     Positive n.
// Returns:
//     Sorted unique primes; 12 -> {2,3}.
// Complexity:
//     O(sqrt n).
// Assumptions:
//     n>=1.
// Notes:
//     Commonly used for inclusion-exclusion, phi and radical.
vector<ll> distinct_prime_factors(ll n) {
	vector<ll> answer; for (auto [prime,exponent] : factorize(n)) {
		(void)exponent; answer.push_back(prime);
	}
	return answer;
}

// Purpose: Count distinct prime factors.
// Input: Positive n.
// Returns: omega(n).
// Complexity: Trial-factorization cost.
// Assumptions: n>=1.
// Notes: Use LinearSieve::omega for many bounded queries.
int omega(ll n) { return factorize(n).size(); }

// Purpose: Count prime factors with multiplicity.
// Input: Positive n.
// Returns: big-Omega(n).
// Complexity: Trial-factorization cost.
// Assumptions: n>=1.
// Notes: Use LinearSieve::bigOmega for many bounded queries.
int big_omega(ll n) { return prime_factors(n).size(); }

// Purpose:
//     Obtain the smallest/largest prime divisor without storing repeated factors.
// Input:
//     Signed-64-bit n.
// Returns:
//     Prime divisor, or -1 for n<=1.
// Complexity:
//     O(sqrt n) worst case.
// Assumptions:
//     None beyond signed range.
// Notes:
//     For bounded queries, smallest factor is simply spf[n].
ll smallest_prime_factor(ll n) {
	if (n <= 1) return -1;
	if (n%2 == 0) return 2;
	for (ll prime = 3; prime <= n/prime; prime += 2) if (n%prime == 0) return prime;
	return n;
}

// Purpose: Find the largest prime divisor of n.
// Input: Signed-64-bit n.
// Returns: Largest prime divisor, or -1 for n<=1.
// Complexity: O(sqrt n) worst case.
// Assumptions: n is moderate enough for trial factorization.
// Notes: Use factorize64(n).back().first for a hard 64-bit input.
ll largest_prime_factor(ll n) {
	if (n <= 1) return -1;
	auto factors = factorize(n); return factors.back().first;
}

// Purpose: Test primality by direct trial division.
// Input: Any signed ll value.
// Returns: True exactly for primes.
// Complexity: O(sqrt n).
// Assumptions: Intended for moderate n.
// Notes: Use is_prime64 for arbitrary 64-bit inputs.
bool is_prime_trial(ll n) {
	if (n < 2) return false;
	return smallest_prime_factor(n) == n;
}

// Purpose:
//     Conveniently compute tau(n) after trial factorization.
// Input:
//     Positive n.
// Returns:
//     Number of positive divisors.
// Complexity:
//     O(sqrt n).
// Assumptions:
//     n>=1.
// Notes:
//     If only tau is needed for many bounded n, use SPF/linear-sieve recurrences.
ll divisor_count(ll n) { return divisor_count(factorize(n)); }

// Purpose:
//     Compute exact sigma from an existing factorization or positive n.
// Input:
//     Sorted factors, or n>=1.
// Returns:
//     Sum of positive divisors.
// Complexity:
//     O(sum exponents), plus factorization for the wrapper.
// Assumptions:
//     Result fits ll.
// Notes:
//     Uses the checked existing implementation and asserts on overflow.
ll divisor_sum(const vector<pair<ll,int>>&factors) {
	auto [fits,value] = divisor_sum_checked(factors); assert(fits); return value;
}

// Purpose: Compute sigma(n) after trial-factorizing n.
// Input: Positive n.
// Returns: Exact positive-divisor sum.
// Complexity: O(sqrt n + sum exponents).
// Assumptions: The result fits ll.
// Notes: Pass an existing factorization to avoid repeated work.
ll divisor_sum(ll n) { return divisor_sum(factorize(n)); }

// Purpose:
//     Compute the product of all positive divisors modulo an arbitrary modulus.
// Input:
//     n>=1, its factorization, and modulus>0.
// Returns:
//     Product modulo modulus.
// Complexity:
//     O(log n + log tau(n)).
// Assumptions:
//     factors reconstruct n.
// Notes:
//     Avoids modular division: pair divisors around n; for odd tau, n is square.
ll divisor_product_mod(ll n, const vector<pair<ll,int>>&factors, ll modulus) {
	assert(n >= 1 && modulus > 0); ll count = divisor_count(factors);
	if (!(count&1)) return power_mod(n,count/2,modulus);
	ll root = integer_sqrt(n); assert((i128)root*root == n);
	return power_mod(root,count,modulus);
}

// Purpose:
//     Compute Euler phi directly from an available factorization.
// Input:
//     n>=1 and its distinct-prime factorization.
// Returns:
//     Number of integers in [1,n] coprime to n.
// Complexity:
//     O(omega(n)).
// Assumptions:
//     factors reconstruct n.
// Notes:
//     Updates result as result-=result/prime.
ll phi_from_factorization(ll n, const vector<pair<ll,int>>&factors) {
	return euler_phi_from_factors(n,factors);
}

// Purpose: Compute Euler's phi for one moderate integer.
// Input: Positive n.
// Returns: Count of integers in [1,n] coprime to n.
// Complexity: O(sqrt n) factorization plus O(omega(n)).
// Assumptions: n>=1.
// Notes: Use a sieve for many bounded phi queries.
ll phi(ll n) { return phi_from_factorization(n,factorize(n)); }

// Purpose:
//     Derive the Mobius function from prime exponents.
// Input:
//     Factorization of a positive integer.
// Returns:
//     0 if any exponent>=2, otherwise (-1)^omega.
// Complexity:
//     O(omega(n)).
// Assumptions:
//     Empty factorization represents 1 and returns 1.
// Notes:
//     This is preferable to a full mu sieve for one/few values.
int mobius_from_factorization(const vector<pair<ll,int>>&factors) {
	for (auto [prime,exponent] : factors) { (void)prime; if (exponent >= 2) return 0; }
	return factors.size()&1 ? -1 : 1;
}

// Purpose: Compute the Mobius function for one moderate integer.
// Input: Positive n.
// Returns: -1, 0, or 1.
// Complexity: O(sqrt n) factorization cost.
// Assumptions: n>=1; mobius(1)=1.
// Notes: Use LinearSieve::mu for a full bounded table.
int mobius(ll n) { return mobius_from_factorization(factorize(n)); }

// Purpose:
//     Compute rad(n), the product of distinct prime divisors.
// Input:
//     Factorization or positive n.
// Returns:
//     Radical; rad(1)=1.
// Complexity:
//     O(omega(n)), plus factorization for wrapper.
// Assumptions:
//     Product fits ll because it divides n.
// Notes:
//     Useful as a prime-support signature.
ll radical(const vector<pair<ll,int>>&factors) {
	ll answer = 1; for (auto [prime,exponent] : factors) { (void)exponent; answer *= prime; }
	return answer;
}

// Purpose: Compute the product of distinct prime divisors of n.
// Input: Positive n.
// Returns: rad(n), with rad(1)=1.
// Complexity: O(sqrt n) factorization cost.
// Assumptions: n>=1.
// Notes: The result always divides n and therefore fits ll.
ll radical(ll n) { return radical(factorize(n)); }

// Purpose:
//     Keep each prime exponent modulo k as a collision-free signature.
// Input:
//     Positive n and k>=1.
// Returns:
//     Sorted {prime,exponent%k} pairs with zero residues omitted.
// Complexity:
//     Factorization cost.
// Assumptions:
//     n>=1.
// Notes:
//     Vectors compare lexicographically and may be used as map keys.
vector<pair<ll,int>> exponent_mod_k(ll n, int k) {
	assert(k >= 1); vector<pair<ll,int>> signature;
	for (auto [prime,exponent] : factorize(n)) if (exponent%k) signature.push_back({prime,exponent%k});
	return signature;
}

// Purpose: Build the missing exponent signature needed to form a kth power.
// Input: A normalized exponent-mod-k signature and k>=1.
// Returns: Sorted complementary prime exponents.
// Complexity: O(signature size).
// Assumptions: Every stored exponent is interpreted modulo k.
// Notes: Multiplying signature and complement gives exponents divisible by k.
vector<pair<ll,int>> kth_power_complement(const vector<pair<ll,int>>&signature, int k) {
	assert(k >= 1); vector<pair<ll,int>> answer;
	for (auto [prime,exponent] : signature) {
		int need = (k-exponent%k)%k; if (need) answer.push_back({prime,need});
	}
	return answer;
}

// Purpose: Represent square-free parity without multiplying primes.
// Input: Sorted factorization.
// Returns: Primes whose exponents are odd.
// Complexity: O(omega(n)).
// Assumptions: Exponents are positive.
// Notes: Safe as a map key even when a numeric kernel would overflow elsewhere.
vector<ll> squarefree_signature(const vector<pair<ll,int>>&factors) {
	vector<ll> answer;
	for (auto [prime,exponent] : factors) if (exponent&1) answer.push_back(prime);
	return answer;
}

// Purpose: Compute the numeric square-free kernel from factor exponents.
// Input: Factorization of a positive ll integer.
// Returns: Product of primes having odd exponent.
// Complexity: O(omega(n)).
// Assumptions: Factors reconstruct an ll value.
// Notes: The kernel divides that value and fits ll.
ll squarefree_kernel(const vector<pair<ll,int>>&factors) {
	i128 answer = 1;
	for (auto [prime,exponent] : factors) if (exponent&1) answer *= prime;
	assert(answer <= LLONG_MAX); return (ll)answer;
}

// Purpose: Recognize a perfect square from factor exponents.
// Input: Sorted factorization.
// Returns: True iff every exponent is even.
// Complexity: O(omega(n)).
// Assumptions: Empty factorization represents 1.
// Notes: Avoids roots and floating point.
bool is_perfect_square_factor(const vector<pair<ll,int>>&factors) {
	for (auto [prime,exponent] : factors) { (void)prime; if (exponent&1) return false; }
	return true;
}

// Purpose: Recognize a moderate positive perfect square through factorization.
// Input: Positive n.
// Returns: True iff n is a square.
// Complexity: O(sqrt n) factorization cost.
// Assumptions: n>=1.
// Notes: integer_sqrt is faster when no factorization is otherwise needed.
bool is_perfect_square_factor(ll n) { return is_perfect_square_factor(factorize(n)); }

// Purpose: Recognize a perfect kth power from factor exponents.
// Input: Factorization and k>=1.
// Returns: True iff every exponent is divisible by k.
// Complexity: O(omega(n)).
// Assumptions: Empty factorization represents 1.
// Notes: Works without reconstructing the integer.
bool is_perfect_kth_power_from_factorization(
		const vector<pair<ll,int>>&factors, int k) {
	assert(k >= 1);
	for (auto [prime,exponent] : factors) { (void)prime; if (exponent%k) return false; }
	return true;
}

// Purpose: Find the largest exponent k for which n=x^k.
// Input: Nonzero signed n other than LLONG_MIN.
// Returns: GCD of exponents, restricted to odd k for negative n.
// Complexity: Trial-factorization cost.
// Assumptions: The convention for n=1 or -1 is 1.
// Notes: Negative perfect powers require an odd exponent.
int max_power_exponent(ll n) {
	if (n == 1 || n == -1) return 1;
	assert(n != LLONG_MIN && n != 0); bool negative = n < 0; if (negative) n = -n;
	int answer = 0;
	for (auto [prime,exponent] : factorize(n)) { (void)prime; answer = nt_gcd(answer,exponent); }
	if (negative) while (!(answer&1)) answer /= 2;
	return answer;
}

// Purpose:
//     Generate all positive divisors from prime exponents.
// Input:
//     Factorization representing a valid positive ll integer.
// Returns:
//     Sorted, duplicate-free divisors.
// Complexity:
//     O(tau(n) log tau(n)) including sorting.
// Assumptions:
//     Every generated divisor fits ll.
// Notes:
//     Sorting is simpler and fast because tau(n) is normally small.
vector<ll> divisors(const vector<pair<ll,int>>&factors) {
	return divisors_from_factors(factors);
}

// Purpose: Factor and enumerate all positive divisors of n.
// Input: Positive moderate n.
// Returns: Sorted duplicate-free divisors.
// Complexity: O(sqrt n + tau(n) log tau(n)).
// Assumptions: n>=1.
// Notes: Pass factors directly when already available.
vector<ll> divisors(ll n) { return divisors(factorize(n)); }

// Purpose: Enumerate all positive divisors except n itself.
// Input: Positive n.
// Returns: Sorted proper divisors; 1 gives an empty vector.
// Complexity: Factorization and divisor-generation cost.
// Assumptions: n>=1.
// Notes: For n>1, divisor 1 is included.
vector<ll> proper_divisors(ll n) {
	auto answer = divisors(n); if (!answer.empty()) answer.pop_back(); return answer;
}

// Purpose: Enumerate divisors shared by two integers.
// Input: Integers whose nonnegative gcd is positive.
// Returns: Sorted divisors of gcd(a,b).
// Complexity: Factorization/divisor-generation cost for the gcd.
// Assumptions: gcd(a,b)>=1.
// Notes: Never factor a and b separately for this query.
vector<ll> common_divisors(ll a, ll b) {
	ll g = nt_gcd(a,b); assert(g >= 1); return divisors(g);
}

// Purpose: Count divisors shared by two integers without enumerating them.
// Input: Integers whose nonnegative gcd is positive.
// Returns: tau(gcd(a,b)).
// Complexity: Trial-factorization cost for the gcd.
// Assumptions: gcd(a,b)>=1.
// Notes: Prefer this over common_divisors when only the count is needed.
ll common_divisor_count(ll a, ll b) {
	ll g = nt_gcd(a,b); assert(g >= 1); return divisor_count(g);
}

// Purpose: Compute the nonnegative gcd of an entire vector.
// Input: Vector of signed ll values.
// Returns: Global gcd; an empty vector gives 0.
// Complexity: O(size*log value).
// Assumptions: Absolute values exclude the LLONG_MIN corner case.
// Notes: Shared divisors of all entries are divisors of this result.
ll gcd_all(const vector<ll>&values) { return nt_gcd(values); }

// Purpose: Enumerate divisors shared by every vector element.
// Input: Vector having positive global gcd.
// Returns: Sorted divisors of gcd_all(values).
// Complexity: GCD plus factorization/divisor generation of the result.
// Assumptions: gcd_all(values)>=1.
// Notes: Factors only one number regardless of vector length.
vector<ll> common_divisors(const vector<ll>&values) {
	ll g = gcd_all(values); assert(g >= 1); return divisors(g);
}

// Purpose: Cache the sorted divisor list for repeated bound/range queries.
// Input: Positive n or its factorization.
// Returns: Methods return counts or nearest divisors; -1 means no neighbor.
// Complexity: Build once; every query is O(log tau(n)).
// Assumptions: Generated divisors fit ll.
// Notes: Use free wrappers for a single query.
struct DivisorList {
	vector<ll> d;
	DivisorList(ll n = 1) : d(divisors(n)) {}
	DivisorList(const vector<pair<ll,int>>&factors) : d(divisors(factors)) {}
	int count_leq(ll x) const { return upper_bound(d.begin(),d.end(),x)-d.begin(); }
	int count_range(ll left, ll right) const {
		return left > right ? 0 : upper_bound(d.begin(),d.end(),right)-lower_bound(d.begin(),d.end(),left);
	}
	ll largest_leq(ll x) const {
		auto it = upper_bound(d.begin(),d.end(),x); return it == d.begin() ? -1 : *prev(it);
	}
	ll smallest_geq(ll x) const {
		auto it = lower_bound(d.begin(),d.end(),x); return it == d.end() ? -1 : *it;
	}
};

// Purpose: Count divisors in a closed interval.
// Input: Positive n and arbitrary left,right bounds.
// Returns: Number of divisors d with left<=d<=right.
// Complexity: One divisor-list build plus O(log tau(n)).
// Assumptions: n>=1.
// Notes: Reuse DivisorList when querying one n many times.
ll count_divisors_in_range(ll n, ll left, ll right) {
	return DivisorList(n).count_range(left,right);
}

// Purpose: Find the greatest divisor not exceeding x.
// Input: Positive n and bound x.
// Returns: Greatest qualifying divisor, or -1.
// Complexity: Divisor generation plus binary search.
// Assumptions: n>=1.
// Notes: Reuse DivisorList for repeated queries.
ll largest_divisor_leq(ll n, ll x) { return DivisorList(n).largest_leq(x); }

// Purpose: Find the least divisor not smaller than x.
// Input: Positive n and bound x.
// Returns: Least qualifying divisor, or -1.
// Complexity: Divisor generation plus binary search.
// Assumptions: n>=1.
// Notes: Reuse DivisorList for repeated queries.
ll smallest_divisor_geq(ll n, ll x) { return DivisorList(n).smallest_geq(x); }

// Purpose:
//     Enumerate unordered factor pairs d*(n/d)=n.
// Input:
//     Positive n.
// Returns:
//     Pairs with first<=second, sorted by first.
// Complexity:
//     O(sqrt n).
// Assumptions:
//     n>=1.
// Notes:
//     Direct scanning is simpler than factorization for one product equation.
vector<pair<ll,ll>> divisor_pairs(ll n) {
	assert(n >= 1); vector<pair<ll,ll>> answer;
	for (ll divisor = 1; divisor <= n/divisor; divisor++)
		if (n%divisor == 0) answer.push_back({divisor,n/divisor});
	return answer;
}

// Purpose:
//     Add/remove a factorization from an exponent accumulator.
// Input:
//     Exponent map, factors and signed multiplicity.
// Returns:
//     Nothing; cnt[p]+=multiplicity*exponent.
// Complexity:
//     O(omega(n) log numberOfStoredPrimes).
// Assumptions:
//     Negative multiplicity is used only when final exponents remain valid.
// Notes:
//     Useful for product/ratio/binomial and sliding multiplicative expressions.
void add_factors(map<ll,ll>&count, const vector<pair<ll,int>>&factors,
		ll multiplicity = 1) {
	for (auto [prime,exponent] : factors) {
		count[prime] += multiplicity*exponent;
		if (!count[prime]) count.erase(prime);
	}
}

// Purpose: Reconstruct an integer from its prime-exponent representation.
// Input: Sorted factors with nonnegative exponents.
// Returns: Success flag and, on success, the reconstructed ll value.
// Complexity: O(sum exponents).
// Assumptions: Primes/exponents describe a nonnegative product.
// Notes: Uses i128 and leaves overflow signaled instead of truncating.
bool reconstruct_factorization(const vector<pair<ll,int>>&factors, ll&value) {
	i128 product = 1;
	for (auto [prime,exponent] : factors) for (int i = 0; i < exponent; i++) {
		product *= prime; if (product > LLONG_MAX) return false;
	}
	value = product; return true;
}

// Purpose: Collapse each p^e factor into one prime-power component.
// Input: Factorization of a positive ll integer.
// Returns: Values p^e in prime order.
// Complexity: O(sum exponents).
// Assumptions: Components fit ll because they divide the represented number.
// Notes: These indivisible components drive coprime factor splits.
vector<ll> prime_powers(const vector<pair<ll,int>>&factors) {
	vector<ll> answer;
	for (auto [prime,exponent] : factors) {
		ll power = 1; for (int i = 0; i < exponent; i++) power *= prime;
		answer.push_back(power);
	}
	return answer;
}

// Purpose:
//     Enumerate x*y=n with gcd(x,y)=1.
// Input:
//     Positive n.
// Returns:
//     All ordered coprime factor pairs.
// Complexity:
//     O(2^omega(n)).
// Assumptions:
//     n fits ll.
// Notes:
//     Each whole prime power is independently assigned to one side.
vector<pair<ll,ll>> coprime_factor_pairs(ll n) {
	assert(n >= 1); auto powers = prime_powers(factorize(n)); vector<pair<ll,ll>> answer;
	for (int mask = 0; mask < (1<<powers.size()); mask++) {
		ll left = 1;
		for (int i = 0; i < (int)powers.size(); i++) if (mask>>i&1) left *= powers[i];
		answer.push_back({left,n/left});
	}
	return answer;
}

// Purpose: Count ordered coprime factor pairs with product n.
// Input: Positive n.
// Returns: 2^omega(n), including both endpoint orders.
// Complexity: Factorization cost.
// Assumptions: n fits positive ll.
// Notes: Every complete prime power independently chooses one endpoint.
ll ordered_coprime_factor_pair_count(ll n) {
	assert(n >= 1); return 1LL<<omega(n);
}

// Purpose: Count unordered coprime factor pairs with product n.
// Input: Positive n.
// Returns: 1 for n=1, otherwise 2^(omega(n)-1).
// Complexity: Factorization cost.
// Assumptions: n fits positive ll.
// Notes: Swapping pair endpoints is considered identical.
ll unordered_coprime_factor_pair_count(ll n) {
	assert(n >= 1); return n == 1 ? 1 : 1LL<<(omega(n)-1);
}

// Purpose:
//     Count divisors d of n that are multiples of x.
// Input:
//     Positive n,x.
// Returns:
//     Product(e_n-e_x+1), or 0 when x does not divide n.
// Complexity:
//     Factorization cost plus O(omega(n)+omega(x)).
// Assumptions:
//     n,x>=1.
// Notes:
//     Each divisor exponent must lie between the exponent in x and that in n.
ll divisors_multiple_of(ll n, ll x) {
	assert(n >= 1 && x >= 1); if (n%x) return 0;
	auto nf = factorize(n), xf = factorize(x); int j = 0; ll answer = 1;
	for (auto [prime,exponent] : nf) {
		int need = j < (int)xf.size() && xf[j].first == prime ? xf[j++].second : 0;
		answer *= exponent-need+1;
	}
	return answer;
}

// Purpose: Count divisors that are perfect kth powers.
// Input: Factorization and k>=1.
// Returns: Product over (floor(exponent/k)+1).
// Complexity: O(omega(n)).
// Assumptions: The count fits ll.
// Notes: Each divisor exponent is chosen from 0,k,2k,... independently.
ll kth_power_divisor_count(const vector<pair<ll,int>>&factors, int k) {
	assert(k >= 1); ll answer = 1;
	for (auto [prime,exponent] : factors) { (void)prime; answer *= exponent/k+1; }
	return answer;
}

// Purpose: Count square divisors.
// Input: Factorization of a positive integer.
// Returns: Product over (floor(exponent/2)+1).
// Complexity: O(omega(n)).
// Assumptions: The count fits ll.
// Notes: This is kth_power_divisor_count with k=2.
ll square_divisor_count(const vector<pair<ll,int>>&factors) {
	return kth_power_divisor_count(factors,2);
}

// Purpose:
//     Construct the largest kth-power divisor of n.
// Input:
//     Positive n and k>=1.
// Returns:
//     Largest divisor whose prime exponents are multiples of k.
// Complexity:
//     Factorization cost.
// Assumptions:
//     Result fits ll because it divides n.
// Notes:
//     For k=2 this is the largest square divisor.
ll largest_kth_power_divisor(ll n, int k) {
	assert(n >= 1 && k >= 1); ll answer = 1;
	for (auto [prime,exponent] : factorize(n))
		for (int i = 0; i < exponent/k*k; i++) answer *= prime;
	return answer;
}

// Purpose: Construct the largest square divisor of n.
// Input: Positive n.
// Returns: Largest divisor whose factor exponents are all even.
// Complexity: Trial-factorization cost.
// Assumptions: n>=1.
// Notes: The result divides n and cannot overflow ll.
ll largest_square_divisor(ll n) { return largest_kth_power_divisor(n,2); }

// Purpose:
//     Find the smallest positive multiplier making n a perfect kth power.
// Input:
//     Positive n and k>=1.
// Returns:
//     Multiplier, or -1 if it does not fit ll.
// Complexity:
//     Factorization cost.
// Assumptions:
//     n>=1.
// Notes:
//     Adds (k-e%k)%k copies of every prime.
ll min_multiplier_to_kth_power(ll n, int k) {
	assert(n >= 1 && k >= 1); i128 answer = 1;
	for (auto [prime,exponent] : factorize(n)) {
		int need = (k-exponent%k)%k;
		for (int i = 0; i < need; i++) { answer *= prime; if (answer > LLONG_MAX) return -1; }
	}
	return answer;
}

// Purpose: Find the smallest multiplier making n a square.
// Input: Positive n.
// Returns: Its square-free kernel.
// Complexity: Trial-factorization cost.
// Assumptions: n>=1.
// Notes: Equal to the smallest divisor to remove so the quotient is square.
ll min_multiplier_to_square(ll n) { return min_multiplier_to_kth_power(n,2); }

// Purpose:
//     Return v_p(n), the number of repeated divisions by p.
// Input:
//     n>=1 and p>=2.
// Returns:
//     Largest exponent e with p^e dividing n.
// Complexity:
//     O(log_p n).
// Assumptions:
//     Usually p is prime, though repeated division is defined for any p>1.
// Notes:
//     For factorial valuations use vp_factorial instead.
int vp(ll n, ll p) {
	assert(n >= 1 && p >= 2); int exponent = 0;
	while (n%p == 0) n /= p, exponent++;
	return exponent;
}

// Purpose: Remove every factor p from n while recording its valuation.
// Input: n>=1 and p>=2.
// Returns: {remaining value,v_p(n)}.
// Complexity: O(log_p n) divisions.
// Assumptions: p is normally prime, though the loop is defined for p>1.
// Notes: The first component is no longer divisible by p.
pair<ll,int> strip_factor(ll n, ll p) {
	assert(n >= 1 && p >= 2); int exponent = 0;
	while (n%p == 0) n /= p, exponent++;
	return {n,exponent};
}

// Purpose:
//     Prime-factorize n! by Legendre valuations.
// Input:
//     n>=0 and a prime list containing every prime<=n.
// Returns:
//     Sorted {prime,exponent} pairs.
// Complexity:
//     O(pi(n) log n).
// Assumptions:
//     Exponents fit ll.
// Notes:
//     Useful when a huge factorial itself must never be constructed.
vector<pair<int,ll>> factorial_factorization(int n, const vector<int>&primes) {
	assert(n >= 0); vector<pair<int,ll>> answer;
	for (int prime : primes) {
		if (prime > n) break;
		answer.push_back({prime,vp_factorial(n,prime)});
	}
	return answer;
}

// Purpose: Prime-factorize C(n,k) using factorial valuations.
// Input: n>=0, integer k, and primes through n.
// Returns: Sorted nonzero {prime,exponent} pairs; invalid k gives empty.
// Complexity: O(pi(n) log n).
// Assumptions: Exponents fit ll.
// Notes: Never constructs factorials or the binomial value.
vector<pair<int,ll>> binomial_factorization(int n, int k,
		const vector<int>&primes) {
	if (k < 0 || k > n) return {};
	vector<pair<int,ll>> answer;
	for (int prime : primes) {
		if (prime > n) break;
		ll exponent = vp_binomial(n,k,prime); if (exponent) answer.push_back({prime,exponent});
	}
	return answer;
}

// Purpose:
//     Compute gcd/lcm/product in prime-exponent representation.
// Input:
//     Two sorted factorization vectors.
// Returns:
//     Sorted factorization using min/max/sum exponents respectively.
// Complexity:
//     O(size(left)+size(right)).
// Assumptions:
//     Inputs contain positive exponents and unique sorted primes.
// Notes:
//     Representation avoids overflow of the underlying integer.
vector<pair<ll,int>> factor_gcd(const vector<pair<ll,int>>&left,
		const vector<pair<ll,int>>&right) {
	vector<pair<ll,int>> answer; int i = 0, j = 0;
	while (i < (int)left.size() && j < (int)right.size()) {
		if (left[i].first < right[j].first) i++;
		else if (right[j].first < left[i].first) j++;
		else answer.push_back({left[i].first,min(left[i].second,right[j].second)}), i++, j++;
	}
	return answer;
}

// Purpose: Compute an lcm in factor-exponent representation.
// Input: Two normalized sorted factorization vectors.
// Returns: Sorted factors with componentwise maximum exponents.
// Complexity: O(size(left)+size(right)).
// Assumptions: Inputs contain unique primes and positive exponents.
// Notes: Does not construct the possibly overflowing lcm value.
vector<pair<ll,int>> factor_lcm(const vector<pair<ll,int>>&left,
		const vector<pair<ll,int>>&right) {
	vector<pair<ll,int>> answer; int i = 0, j = 0;
	while (i < (int)left.size() || j < (int)right.size()) {
		if (j == (int)right.size() || (i < (int)left.size() && left[i].first < right[j].first)) answer.push_back(left[i++]);
		else if (i == (int)left.size() || right[j].first < left[i].first) answer.push_back(right[j++]);
		else answer.push_back({left[i].first,max(left[i].second,right[j].second)}), i++, j++;
	}
	return answer;
}

// Purpose: Multiply two factorized integers symbolically.
// Input: Two normalized sorted factorization vectors.
// Returns: Sorted factors with componentwise summed exponents.
// Complexity: O((L+R) log(L+R)) in this compact implementation.
// Assumptions: Summed exponents fit int.
// Notes: The represented numeric product need not fit ll.
vector<pair<ll,int>> factor_multiply(const vector<pair<ll,int>>&left,
		const vector<pair<ll,int>>&right) {
	vector<pair<ll,int>> answer = factor_lcm(left,right);
	for (auto&[prime,exponent] : answer) {
		auto a = lower_bound(left.begin(),left.end(),pair<ll,int>{prime,-1});
		auto b = lower_bound(right.begin(),right.end(),pair<ll,int>{prime,-1});
		exponent = (a != left.end() && a->first == prime ? a->second : 0)
			+(b != right.end() && b->first == prime ? b->second : 0);
	}
	return answer;
}

// Purpose: Divide one factorized integer by another symbolically.
// Input: Numerator and denominator factorizations.
// Returns: Sorted difference of exponent vectors.
// Complexity: O((L+R) log(L+R)).
// Assumptions: Denominator divides numerator and result exponents fit int.
// Notes: Assertions reject negative final exponents.
vector<pair<ll,int>> factor_divide(const vector<pair<ll,int>>&numerator,
		const vector<pair<ll,int>>&denominator) {
	map<ll,ll> exponents; add_factors(exponents,numerator); add_factors(exponents,denominator,-1);
	vector<pair<ll,int>> answer;
	for (auto [prime,exponent] : exponents) { assert(exponent >= 0); if (exponent) answer.push_back({prime,(int)exponent}); }
	return answer;
}

// Purpose: Test divisibility using only prime-exponent vectors.
// Input: Candidate divisor and number factorizations.
// Returns: True iff every divisor exponent is available in the number.
// Complexity: O(size(divisor)+size(number)).
// Assumptions: Both inputs are normalized and sorted.
// Notes: Avoids reconstructing either integer.
bool factor_divides(const vector<pair<ll,int>>&divisor,
		const vector<pair<ll,int>>&number) {
	int j = 0;
	for (auto [prime,exponent] : divisor) {
		while (j < (int)number.size() && number[j].first < prime) j++;
		if (j == (int)number.size() || number[j].first != prime || number[j].second < exponent) return false;
	}
	return true;
}

// Purpose: Compute a divisor count modulo an arbitrary modulus.
// Input: Factorization and modulus>0.
// Returns: Product(exponent+1) modulo modulus.
// Complexity: O(omega(n)).
// Assumptions: mul_mod supports the supplied positive modulus.
// Notes: Useful when the represented product or its tau value is huge.
ll divisor_count_mod(const vector<pair<ll,int>>&factors, ll modulus) {
	assert(modulus > 0); ll answer = 1%modulus;
	for (auto [prime,exponent] : factors) { (void)prime; answer = mul_mod(answer,exponent+1,modulus); }
	return answer;
}

// Purpose:
//     Enumerate square-free divisors and their Mobius signs for one fixed n.
// Input:
//     Factorization of n.
// Returns:
//     squarefree_divisors: values; mobius_divisors: {value,mu(value)}.
// Complexity:
//     O(2^omega(n) log 2^omega(n)) including sorting.
// Assumptions:
//     Products divide n and fit ll.
// Notes:
//     Replaces a global Mobius sieve when only divisors of one huge n matter.
vector<pair<ll,int>> mobius_divisors(const vector<pair<ll,int>>&factors) {
	vector<pair<ll,int>> answer{{1,1}};
	for (auto [prime,exponent] : factors) {
		(void)exponent; int old = answer.size();
		for (int i = 0; i < old; i++) answer.push_back({answer[i].first*prime,-answer[i].second});
	}
	sort(answer.begin(),answer.end()); return answer;
}

// Purpose: Enumerate square-free divisors of one factorized integer.
// Input: Factorization of n.
// Returns: Sorted divisors formed by choosing every distinct prime 0/1 times.
// Complexity: O(2^omega(n) log 2^omega(n)).
// Assumptions: Products fit ll because they divide n.
// Notes: Useful for divisor inclusion-exclusion.
vector<ll> squarefree_divisors(const vector<pair<ll,int>>&factors) {
	auto withMu = mobius_divisors(factors); vector<ll> answer;
	for (auto [value,mu] : withMu) { (void)mu; answer.push_back(value); }
	return answer;
}

// Purpose: Enumerate safe subset products with inclusion-exclusion signs.
// Input: Distinct primes and a positive product limit.
// Returns: {product,(-1)^subsetSize}; products above the limit are pruned.
// Complexity: O(number of valid subsets), at most O(2^k).
// Assumptions: At most 25 input primes.
// Notes: Product 1 has sign +1; overflow is prevented by division checks.
vector<pair<ll,int>> subset_products_with_parity(const vector<ll>&primes,
		ll limit) {
	assert(limit >= 1 && primes.size() <= 25); vector<pair<ll,int>> answer{{1,1}};
	for (ll prime : primes) {
		int old = answer.size();
		for (int i = 0; i < old; i++) if (answer[i].first <= limit/prime)
			answer.push_back({answer[i].first*prime,-answer[i].second});
	}
	return answer;
}

// Purpose:
//     Count values in prefixes/ranges with a specified gcd against fixed n.
// Input:
//     Bounds, positive n, and optional gcd target g.
// Returns:
//     Coprime or exact-gcd counts.
// Complexity:
//     O(2^omega(n)) after trial factorization.
// Assumptions:
//     Bounds are interpreted over positive x; n>=1.
// Notes:
//     gcd(x,n)=g transforms to x=g*y and gcd(y,n/g)=1.
ll count_coprime_leq(ll maximum, ll n) {
	return maximum <= 0 ? 0 : count_coprime_to_fixed(maximum,n);
}

// Purpose: Count positive integers in a closed range coprime to fixed n.
// Input: Bounds left,right and n>=1.
// Returns: Count of x in [left,right] with gcd(x,n)=1.
// Complexity: O(2^omega(n)) per prefix, plus factorization.
// Assumptions: Only positive x are counted.
// Notes: Empty/reversed intervals return 0.
ll count_coprime_range(ll left, ll right, ll n) {
	return left > right ? 0 : count_coprime_leq(right,n)-count_coprime_leq(left-1,n);
}

// Purpose: Count positive x<=maximum with gcd(x,n) equal to a target.
// Input: Prefix bound, n>=1 and gcdValue>=1.
// Returns: Exact-gcd count, or 0 if gcdValue does not divide n.
// Complexity: O(2^omega(n/gcdValue)) plus factorization.
// Assumptions: Only positive x are counted.
// Notes: Divide x and n by gcdValue, then impose coprimality.
ll count_gcd_eq_leq(ll maximum, ll n, ll gcdValue) {
	assert(n >= 1 && gcdValue >= 1); if (n%gcdValue) return 0;
	return count_coprime_leq(maximum/gcdValue,n/gcdValue);
}

// Purpose: Count divisors of n that are coprime to another integer.
// Input: n>=1 and arbitrary other.
// Returns: Product(exponent+1) over primes of n absent from other.
// Complexity: Trial-factorization cost.
// Assumptions: n>=1.
// Notes: A forbidden prime forces its divisor exponent to zero.
ll count_divisors_coprime_to(ll n, ll other) {
	assert(n >= 1); ll answer = 1;
	for (auto [prime,exponent] : factorize(n)) if (other%prime) answer *= exponent+1;
	return answer;
}

// Purpose: Sum divisors of n that are coprime to another integer.
// Input: n>=1 and arbitrary other.
// Returns: Product of allowed geometric prime-power sums.
// Complexity: Trial-factorization plus O(sum exponents).
// Assumptions: The exact result fits ll.
// Notes: Prime powers shared with other contribute only exponent zero.
ll sum_divisors_coprime_to(ll n, ll other) {
	assert(n >= 1); i128 answer = 1;
	for (auto [prime,exponent] : factorize(n)) if (other%prime) {
		i128 term = 1, power = 1;
		for (int i = 0; i < exponent; i++) power *= prime, term += power;
		answer *= term;
	}
	assert(answer <= LLONG_MAX); return answer;
}

// Purpose:
//     Sum a bounded frequency/value array over multiples of every divisor.
// Input:
//     values[0..M], where index zero is ignored.
// Returns:
//     result[d]=sum(values[m]) over positive multiples m of d.
// Complexity:
//     O(M log M).
// Assumptions:
//     Additions fit ll.
// Notes:
//     With a frequency array this counts input numbers divisible by d.
vector<ll> multiple_zeta(const vector<ll>&values) {
	int maximum = values.size()-1; vector<ll> answer(maximum+1);
	for (int divisor = 1; divisor <= maximum; divisor++)
		for (int multiple = divisor; multiple <= maximum; multiple += divisor) answer[divisor] += values[multiple];
	return answer;
}

// Purpose: Give the frequency-array operation its common contest name.
// Input: frequency[0..M], index zero ignored.
// Returns: result[d]=number/sum of values at multiples of d.
// Complexity: O(M log M).
// Assumptions: Sums fit ll.
// Notes: This is an alias of multiple_zeta.
vector<ll> count_multiples(const vector<ll>&frequency) {
	return multiple_zeta(frequency);
}

// Purpose: Aggregate an indexed function over divisors.
// Input: values[0..M], index zero ignored.
// Returns: result[x]=sum(values[d]) for every d dividing x.
// Complexity: O(M log M).
// Assumptions: Sums fit ll.
// Notes: Direction is divisors into each multiple; inverse uses divisor_mobius_transform.
vector<ll> divisor_zeta(vector<ll> values) {
	return divisor_zeta_transform(move(values));
}

// Purpose:
//     Generate only divisors no larger than a supplied bound.
// Input:
//     Factorization and nonnegative limit.
// Returns:
//     Sorted qualifying divisors.
// Complexity:
//     O(number of explored valid partial products + output log output).
// Assumptions:
//     Prime factors are positive and sorted.
// Notes:
//     Checks current<=limit/prime before multiplying.
vector<ll> divisors_leq(const vector<pair<ll,int>>&factors, ll limit) {
	if (limit < 1) return {};
	vector<ll> answer;
	function<void(int,ll)> generate = [&](int at,ll current) {
		if (at == (int)factors.size()) { answer.push_back(current); return; }
		auto [prime,exponent] = factors[at]; ll value = current;
		for (int e = 0; e <= exponent; e++) {
			generate(at+1,value); if (e == exponent || value > limit/prime) break; value *= prime;
		}
	};
	generate(0,1); sort(answer.begin(),answer.end()); return answer;
}

// Purpose: Bundle bounded SPF factorization with phi, mu and divisor queries.
// Input: Maximum precomputation bound.
// Returns: Member methods provide factor data for values inside that bound.
// Complexity: O(N) build; factor queries O(bigOmega(x)); divisor output O(tau(x)).
// Assumptions: Values queried satisfy 1<=x<=maximum.
// Notes: Prefer this mode for many bounded values; it intentionally reuses LinearSieve.
struct FactorSieve {
	LinearSieve data;
	FactorSieve(int maximum) : data(maximum) {}

	vector<pair<int,int>> factorize(int value) const {
		auto factors = factorize_spf(value,data.spf); vector<pair<int,int>> answer;
		for (auto [prime,exponent] : factors) answer.push_back({(int)prime,exponent});
		return answer;
	}
	vector<int> distinct_factors(int value) const {
		vector<int> answer; for (auto [prime,exponent] : factorize(value)) {
			(void)exponent; answer.push_back(prime);
		}
		return answer;
	}
	vector<int> divisors(int value) const {
		vector<pair<ll,int>> factors;
		for (auto [prime,exponent] : factorize(value)) factors.push_back({prime,exponent});
		auto values = divisors_from_factors(factors); return vector<int>(values.begin(),values.end());
	}
	int phi(int value) const { assert(1 <= value && value <= data.n); return data.phi[value]; }
	int mobius(int value) const { assert(1 <= value && value <= data.n); return data.mu[value]; }
};

/*
DIVISOR / FACTORIZATION HELPER CATALOG

Choose a factorizer:
  factorize_spf(n,spf)          many n inside a precomputed bound
  factorize(n)                  one/few moderate n, no precomputation
  factorize(n,primes)           one/few n with primes through sqrt(n)
  factorize64(n)                arbitrary positive signed-64-bit n

Factor representations:
  factorize*, FactorSieve       sorted {prime,exponent} pairs
  prime_factors*                primes repeated by multiplicity
  distinct_prime_factors       every prime divisor once
  omega / big_omega             distinct / repeated factor counts

Divisors and multiplicative functions:
  divisors, proper_divisors, divisor_pairs, DivisorList
  divisor_count/sum/product_mod, phi, mobius, radical
  squarefree_kernel, exponent_mod_k, kth_power_complement
  kth_power_divisor_count, largest_kth_power_divisor

Factor algebra and valuations:
  factor_gcd/lcm/multiply/divide/divides, reconstruct_factorization
  vp, strip_factor, factorial_factorization, binomial_factorization

Counting and transforms:
  count_coprime_leq/range, count_gcd_eq_leq
  count_divisors_coprime_to, sum_divisors_coprime_to
  divisor_zeta, multiple_zeta/count_multiples

Important constraints:
  Passing an incomplete prime list to factorize(n,primes) is invalid.
  Explicit divisor generation costs Theta(tau(n)) memory and time.
  Pollard Rho is randomized and intended for the genuinely 64-bit path.
  Product/sum helpers assert or return -1 where their documented ll bound fails.
*/

ll radical_from_factors(const vector<pair<ll,int>>&factors) {
	return radical(factors);
}

int maximum_perfect_power_exponent(ll value) {
	assert(value != LLONG_MIN && value != 0 && value != 1 && value != -1);
	bool negative = value < 0;
	ll positive = negative ? -value : value; int answer = 0;
	for (auto [prime,exponent] : factorize_64(positive)) {
		(void)prime; answer = nt_gcd(answer,exponent);
	}
	if (negative) while (!(answer&1)) answer /= 2;
	return answer;
}

// Problem:
//     Find the maximum length of a subsequence reorderable as x1|x2|... .
// Mathematical idea:
//     DP from each divisor to its multiples; all duplicate copies can be consecutive.
// State / Precomputation:
//     Frequency by bounded value and best chain ending at each value.
// Returns:
//     Maximum selected element count.
// Complexity:
//     O(M log M+n).
// Assumptions:
//     Positive bounded values; reordering is allowed.
// Notes:
//     This is divisibility-poset DP, not LIS in the original order.
int longest_divisibility_chain(const vector<int>&values) {
	if (values.empty()) return 0;
	int maximum = *max_element(values.begin(),values.end());
	vector<int> frequency(maximum+1),dp(maximum+1); for (int x : values) frequency[x]++;
	int answer = 0;
	for (int divisor = 1; divisor <= maximum; divisor++) {
		dp[divisor] = max(dp[divisor],frequency[divisor]); answer = max(answer,dp[divisor]);
		for (int multiple = divisor+divisor; multiple <= maximum; multiple += divisor)
			dp[multiple] = max(dp[multiple],dp[divisor]+frequency[multiple]);
	}
	return answer;
}


#ifndef NUMBER_THEORY_TEMPLATE_NO_MAIN
int main() {
	ios::sync_with_stdio(false);
	cin.tie(nullptr);

	ll x,y,g;
	if (find_any_solution(15,21,6,x,y,g)) cout << x << ' ' << y << '\n';

	auto [remainder,modulus] = crt_merge(2,6,5,9);
	cout << remainder << " modulo " << modulus << '\n';

	LinearSieve sieve(1000000);
	for (auto [prime,exponent] : factorize_spf(360,sieve.spf))
		cout << prime << '^' << exponent << ' ';
	cout << '\n';

	cout << divisor_count(360) << '\n';
	cout << power_mod(2,100,1000000007) << '\n';
	cout << mod_sqrt_prime(10,13) << '\n';
}
#endif
