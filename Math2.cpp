/*
COMBINATORICS -- SAM-STYLE C++17/C++20 PRINTABLE TEMPLATE

CONTAINS
- A self-contained arithmetic foundation: GCD, modular operations, CRT,
  deterministic 64-bit primality, and Pollard-Rho factorization.
- Prime-mod nCr/nPr with factorial and inverse-factorial tables.
- Lucas theorem for huge n under a small prime modulus.
- Legendre valuations for n!, binomial coefficients, and trailing zeros.
- Binomial coefficients modulo p^q and manageable composite moduli.
- Bounded stars and bars by inclusion-exclusion.
- Prime factorizations of n! and C(n,k).
- Count of nonzero binomial coefficients modulo a prime.

IMPORTANT CHOICES
- Prime modulus and n<mod: Comb.
- Huge n and small prime modulus: lucas.
- Composite modulus with manageable prime powers: ncr_mod_composite.
- Only divisibility by p^e: vp_binomial; do not construct C(n,k).

BUILD
    g++ -std=c++17 -O2 combinatorics_template.cpp

STRUCTURE
    Combinatorics cb;
    Combinatorics::Comb comb(N,mod); // reusable factorial state
    comb.C(n,r);                     // O(1) prime-mod query
    cb.lucas(n,r,p);                 // stateless algorithms are methods
*/

#include <bits/stdc++.h>
using namespace std;

using ll = long long;
__extension__ typedef __int128 i128;
__extension__ typedef unsigned __int128 u128;

struct Combinatorics {
	// Stateless facade; Comb/PrimePowerBinomial own precomputed state. O(1)
	Combinatorics() = default;
	
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
	static ll nt_gcd(ll a, ll b) {
		a = abs(a); b = abs(b);
		while (b) { ll r = a%b; a = b; b = r; }
		return a;
	}
	
	static ll nt_gcd(const vector<ll>&a) {
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
	static pair<bool,ll> nt_lcm_checked(ll a, ll b) {
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
	static ll floor_div(ll a, ll b) {
		assert(b != 0); ll q = a/b, r = a%b;
		if (r && ((r > 0) != (b > 0))) q--;
		return q;
	}
	
	static ll ceil_div(ll a, ll b) {
		assert(b != 0); ll q = a/b, r = a%b;
		if (r && ((r > 0) == (b > 0))) q++;
		return q;
	}
	
	static i128 floor_div128(i128 a, i128 b) {
		assert(b != 0); i128 q = a/b, r = a%b;
		if (r && ((r > 0) != (b > 0))) q--;
		return q;
	}
	
	static i128 ceil_div128(i128 a, i128 b) {
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
	static ll extgcd(ll a, ll b, ll&x, ll&y) {
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
	static bool find_any_solution(ll a, ll b, ll c, ll&x, ll&y, ll&g) {
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
	static ll count_diophantine_solutions(ll a, ll b, ll c,
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
	
	
	static ll mod_normalize(ll x, ll mod) {
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
	static ll mul_mod(ll a, ll b, ll mod) {
		return (i128)mod_normalize(a,mod)*mod_normalize(b,mod)%mod;
	}
	
	static ll power_mod(ll a, unsigned long long exponent, ll mod) {
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
	static ll inverse_mod(ll a, ll mod) {
		assert(mod > 1); ll x,y,g = extgcd(a,mod,x,y);
		return g == 1 ? mod_normalize(x,mod) : -1;
	}
	
	static vector<ll> inverses_prime(int n, ll primeMod) {
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
	static pair<ll,ll> power_geometric_sum(ll a, unsigned long long n, ll mod) {
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
	static pair<ll,ll> crt_merge(ll a, ll m, ll b, ll n) {
		assert(m > 0 && n > 0); a = mod_normalize(a,m); b = mod_normalize(b,n);
		ll x,y,g = extgcd(m,n,x,y); ll difference = b-a;
		if (difference%g) return {0,-1};
		i128 modulus = (i128)(m/g)*n; if (modulus > LLONG_MAX) return {0,-1};
		ll reduced = n/g;
		ll k = mul_mod(difference/g,mod_normalize(x,reduced),reduced);
		ll result = (ll)(((i128)a+(i128)m*k)%modulus);
		return {result,(ll)modulus};
	}
	
	static pair<ll,ll> crt_system(const vector<ll>&remainders, const vector<ll>&moduli) {
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
	static pair<ll,ll> solve_linear_congruence(ll a, ll b, ll modulus) {
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
	static pair<vector<int>,vector<char>> eratosthenes(int n) {
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
	static vector<pair<int,int>> factorize_spf(int x, const vector<int>&spf) {
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
	static vector<pair<ll,int>> factorize_trial(ll n, const vector<int>&primes) {
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
	static bool is_prime_64(ll n) {
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
	
	inline static mt19937_64 pollard_rng{
		(unsigned long long)chrono::steady_clock::now().time_since_epoch().count()};
	
	static ll pollard_rho_factor(ll n) {
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
	
	static void factorize_64_raw(ll n, vector<ll>&factors) {
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
	static vector<pair<ll,int>> factorize_64(ll n) {
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
	// Complexity:
	//     Roughly O(sqrt(R) log log R + (R-L+1) log log R).
	// Assumptions:
	//     0<=left<=right and the segment fits memory.
	// Notes:
	//     Start at max(p^2,ceil(left/p)*p), so p itself is not marked composite.
	static vector<ll> segmented_primes(ll left, ll right) {
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
	
	
	struct Comb {
		int maximum;
		ll mod;
		vector<ll> factorial,inverseFactorial;
	
		// Problem:
		//     Answer nCr and nPr quickly under a prime modulus.
		// Mathematical idea:
		//     nCr=n!/(r!(n-r)!), with Fermat inverses because mod is prime.
		// State / Precomputation:
		//     factorial and inverse factorial through maximum.
		// Returns:
		//     C(n,r) or P(n,r) modulo mod.
		// Complexity:
		//     O(maximum+log mod) build, O(1) query.
		// Assumptions:
		//     mod is prime and maximum<mod.
		// Notes:
		//     If n>=mod factorial becomes zero; use Lucas or p-adic methods instead.
		Comb(int limit = 0, ll primeMod = 1000000007)
			: maximum(limit), mod(primeMod), factorial(limit+1,1), inverseFactorial(limit+1,1) {
			assert(0 <= limit && limit < mod && is_prime_64(mod));
			for (int i = 1; i <= limit; i++) factorial[i] = mul_mod(factorial[i-1],i,mod);
			if (limit) inverseFactorial[limit] = power_mod(factorial[limit],mod-2,mod);
			for (int i = limit; i >= 1; i--) inverseFactorial[i-1] = mul_mod(inverseFactorial[i],i,mod);
		}
	
		ll C(ll n, ll r) const {
			if (r < 0 || n < r || n < 0 || n > maximum) return 0;
			return mul_mod(factorial[n],mul_mod(inverseFactorial[r],inverseFactorial[n-r],mod),mod);
		}
		ll P(ll n, ll r) const {
			if (r < 0 || n < r || n < 0 || n > maximum) return 0;
			return mul_mod(factorial[n],inverseFactorial[n-r],mod);
		}
	};
	
	// Problem:
	//     Compute C(n,r) modulo a small prime when n may be huge.
	// Mathematical idea:
	//     Lucas: multiply C(base-p digit of n,base-p digit of r).
	// State / Precomputation:
	//     Factorials through p-1.
	// Returns:
	//     C(n,r) mod p.
	// Complexity:
	//     O(p+log_p n).
	// Assumptions:
	//     p is prime and small enough for O(p) memory.
	// Notes:
	//     A digit r_i>n_i makes the answer zero.
	static ll lucas(ll n, ll r, int prime) {
		assert(n >= 0 && r >= 0 && is_prime_64(prime)); Comb comb(prime-1,prime); ll answer = 1;
		while (n || r) {
			int ni = n%prime, ri = r%prime; if (ri > ni) return 0;
			answer = mul_mod(answer,comb.C(ni,ri),prime); n /= prime; r /= prime;
		}
		return answer;
	}
	
	static ll vp_factorial(ll n, ll prime) {
		assert(n >= 0 && prime >= 2); ll exponent = 0;
		while (n) n /= prime, exponent += n;
		return exponent;
	}
	
	static ll vp_binomial(ll n, ll r, ll prime) {
		if (r < 0 || r > n) return 0;
		return vp_factorial(n,prime)-vp_factorial(r,prime)-vp_factorial(n-r,prime);
	}
	
	static ll trailing_zeros_factorial_base(ll n, ll base) {
		assert(n >= 0 && base >= 2); ll answer = LLONG_MAX;
		for (auto [prime,exponent] : factorize_64(base)) answer = min(answer,vp_factorial(n,prime)/exponent);
		return answer;
	}
	
	struct PrimePowerBinomial {
		ll prime,mod;
		int exponent;
		vector<ll> prefix;
	
		PrimePowerBinomial(ll primeValue, int exponentValue)
			: prime(primeValue), mod(1), exponent(exponentValue) {
			assert(is_prime_64(prime) && exponent >= 1);
			for (int i = 0; i < exponent; i++) {
				assert((i128)mod*prime <= 3000000); mod *= prime;
			}
			prefix.assign(mod+1,1);
			for (ll i = 1; i <= mod; i++) {
				prefix[i] = prefix[i-1];
				if (i%prime) prefix[i] = mul_mod(prefix[i],i,mod);
			}
		}
	
		ll unit_factorial(ll n) const {
			ll answer = 1;
			while (n > 1) {
				answer = mul_mod(answer,power_mod(prefix[mod],n/mod,mod),mod);
				answer = mul_mod(answer,prefix[n%mod],mod); n /= prime;
			}
			return answer;
		}
	
		// Problem:
		//     Compute C(n,r) modulo p^q without illegal factorial division.
		// Mathematical idea:
		//     Separate p-adic valuation from the factorial unit parts.
		// State / Precomputation:
		//     Products excluding multiples of p over one period p^q.
		// Returns:
		//     C(n,r) mod p^q.
		// Complexity:
		//     O(p^q) build and O(log_p n+log mod) query.
		// Assumptions:
		//     p is prime and p^q<=3e6 in this implementation.
		// Notes:
		//     Unit factorials are coprime to p^q and therefore invertible.
		ll C(ll n, ll r) const {
			if (r < 0 || r > n) return 0;
			ll valuation = vp_binomial(n,r,prime); if (valuation >= exponent) return 0;
			ll numerator = unit_factorial(n), left = unit_factorial(r), right = unit_factorial(n-r);
			ll answer = mul_mod(numerator,inverse_mod(left,mod),mod);
			answer = mul_mod(answer,inverse_mod(right,mod),mod);
			return mul_mod(answer,power_mod(prime,valuation,mod),mod);
		}
	};
	
	// Problem:
	//     Compute C(n,r) modulo a manageable arbitrary composite modulus.
	// Mathematical idea:
	//     Solve modulo each coprime prime power, then combine by CRT.
	// State / Precomputation:
	//     One PrimePowerBinomial per factor p^q.
	// Returns:
	//     C(n,r) mod modulus.
	// Complexity:
	//     O(sum p^q + numberOfFactors*log n).
	// Assumptions:
	//     Every prime-power factor is <=3e6; modulus fits ll.
	// Notes:
	//     For small n and arbitrary modulus, Pascal DP may be simpler.
	static ll ncr_mod_composite(ll n, ll r, ll modulus) {
		assert(modulus >= 1); if (modulus == 1 || r < 0 || r > n) return 0;
		pair<ll,ll> answer{0,1};
		for (auto [prime,exponent] : factorize_64(modulus)) {
			PrimePowerBinomial comb(prime,exponent);
			answer = crt_merge(answer.first,answer.second,comb.C(n,r),comb.mod);
			assert(answer.second != -1);
		}
		return answer.first;
	}
	
	// Problem:
	//     Count solutions x1+...+variables=sum with 0<=xi<=upper.
	// Mathematical idea:
	//     Stars and bars, excluding subsets forced above upper by inclusion-exclusion.
	// State / Precomputation:
	//     Prime-modulus combination table covering required n.
	// Returns:
	//     Count modulo comb.mod.
	// Complexity:
	//     O(variables).
	// Assumptions:
	//     variables>=0,upper>=0 and required C arguments fit comb.maximum.
	// Notes:
	//     Without upper bound the answer is C(sum+variables-1,variables-1).
	static ll bounded_stars_and_bars(ll sum, int variables, ll upper, const Comb&comb) {
		if (variables == 0) return sum == 0;
		ll answer = 0;
		for (int chosen = 0; chosen <= variables; chosen++) {
			ll remaining = sum-(ll)chosen*(upper+1); if (remaining < 0) break;
			ll term = mul_mod(comb.C(variables,chosen),comb.C(remaining+variables-1,variables-1),comb.mod);
			answer = mod_normalize(answer+(chosen&1 ? -term : term),comb.mod);
		}
		return answer;
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
	static vector<pair<int,ll>> factorial_factorization(int n, const vector<int>&primes) {
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
	static vector<pair<int,ll>> binomial_factorization(int n, int k,
			const vector<int>&primes) {
		if (k < 0 || k > n) return {};
		vector<pair<int,ll>> answer;
		for (int prime : primes) {
			if (prime > n) break;
			ll exponent = vp_binomial(n,k,prime); if (exponent) answer.push_back({prime,exponent});
		}
		return answer;
	}
	
	// Problem:
	//     Count k in [0,n] for which C(n,k) is nonzero modulo prime p.
	// Mathematical idea:
	//     Lucas says every base-p digit of k may independently be 0..n_i.
	// State / Precomputation:
	//     Base-p digits of n.
	// Returns:
	//     Product of (n_i+1).
	// Complexity:
	//     O(log_p n).
	// Assumptions:
	//     p is prime and answer fits ll.
	// Notes:
	//     Equivalently, no base-p carry occurs in k+(n-k) (Kummer viewpoint).
	static ll count_nonzero_binomials_mod_prime(ll n, ll prime) {
		assert(n >= 0 && is_prime_64(prime)); i128 answer = 1;
		while (n) answer *= n%prime+1, n /= prime;
		assert(answer <= LLONG_MAX); return answer;
	}

};

#ifndef COMBINATORICS_TEMPLATE_NO_MAIN
int main() {
	ios::sync_with_stdio(false);
	cin.tie(nullptr);
	Combinatorics cb;

	Combinatorics::Comb comb(1000,1000000007);
	cout << comb.C(10,3) << '\n';
	cout << comb.P(10,3) << '\n';
	cout << cb.lucas(1000000000000LL,12345,101) << '\n';
	cout << cb.ncr_mod_composite(20,7,72) << '\n';

	auto [primes,isPrime] = cb.eratosthenes(100);
	(void)isPrime;
	for (auto [prime,exponent] : cb.binomial_factorization(20,7,primes))
		cout << prime << '^' << exponent << ' ';
	cout << '\n';

	cout << cb.count_nonzero_binomials_mod_prime(100,5) << '\n';
}
#endif
