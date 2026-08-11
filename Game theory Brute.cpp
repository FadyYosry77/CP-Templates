#define GAME_THEORY_TEMPLATE_NO_MAIN
#include "game_theory_template.cpp"

/*
GameTheoryLab: brute force first, formula second.

Quick guide:
1. Write moves(state) exactly; terminal states return an empty vector.
2. Keep the brute limit tiny and make every move decrease a clear measure.
3. Print losing states and Grundy values, then annotate parity, residues, XOR,
   gaps, binary digits, or prime exponents.
4. State a candidate formula and use firstCounterexample... to attack it.
5. Increase the limit and test boundary-shaped states, not only random states.
6. Only after it survives, prove: every move leaves the proposed P-set and
   every position outside it has a move into the P-set.

Build the included subtraction demo:
g++ -std=c++17 -O2 game_theory_bruteforce_lab.cpp

This file is intentionally separate from the contest template. It uses maps and
materializes moves for clarity; do not copy the laboratory solver into a tight
submission when a proved O(1) or iterative formula is available.
*/

struct GameTheoryLab {
	template<class State>
	struct Solver {
		function<vector<State>(const State&)> moves;
		map<State,char> winMemo;
		map<State,int> grundyMemo;
		set<State> activeWin, activeGrundy;

		explicit Solver(function<vector<State>(const State&)> generator)
			: moves(move(generator)) {}

		// Normal play: terminal is losing. Throws if a reachable cycle is found.
		bool win(const State&state) {
			auto found = winMemo.find(state);
			if (found != winMemo.end()) return found->second;
			if (!activeWin.insert(state).second)
				throw invalid_argument("GameTheoryLab::Solver::win: cyclic game");
			bool answer = false;
			for (const State&to : moves(state)) if (!win(to)) {
				answer = true; break;
			}
			activeWin.erase(state);
			winMemo[state] = answer;
			return answer;
		}

		// Sprague-Grundy value for a finite acyclic impartial normal-play game.
		int grundy(const State&state) {
			auto found = grundyMemo.find(state);
			if (found != grundyMemo.end()) return found->second;
			if (!activeGrundy.insert(state).second)
				throw invalid_argument("GameTheoryLab::Solver::grundy: cyclic game");
			vector<int> values;
			for (const State&to : moves(state)) values.push_back(grundy(to));
			activeGrundy.erase(state);
			return grundyMemo[state] = GameTheory::mex(values);
		}

		void clear() {
			winMemo.clear(); grundyMemo.clear();
			activeWin.clear(); activeGrundy.clear();
		}
	};

	static vector<int> canonicalPiles(vector<int> piles, bool eraseZero = true) {
		for (int x : piles) assert(x >= 0);
		if (eraseZero) piles.erase(remove(piles.begin(),piles.end(),0),piles.end());
		sort(piles.begin(),piles.end());
		return piles;
	}

	// Moves for one heap and a fixed subtraction set.
	static function<vector<int>(const int&)> subtractionMoves(vector<int> allowed) {
		sort(allowed.begin(),allowed.end());
		allowed.erase(unique(allowed.begin(),allowed.end()),allowed.end());
		for (int d : allowed) assert(d > 0);
		return [allowed](const int&x) {
			vector<int> out;
			for (int d : allowed) if (d <= x) out.push_back(x-d);
			return out;
		};
	}

	// CSES Grundy's Game: split one chosen heap into unequal positive heaps.
	static vector<vector<int>> unequalSplitMoves(const vector<int>&state) {
		vector<vector<int>> out;
		for (int i = 0; i < (int)state.size(); ++i)
			for (int a = 1; a < state[i]-a; ++a) {
				int b = state[i]-a;
				if (a == b) continue;
				vector<int> to = state;
				to.erase(to.begin()+i); to.push_back(a); to.push_back(b);
				out.push_back(canonicalPiles(to));
			}
		sort(out.begin(),out.end()); out.erase(unique(out.begin(),out.end()),out.end());
		return out;
	}

	// Ordinary Nim moves for small states. maxNewSize can cap generated targets.
	static vector<vector<int>> nimMoves(const vector<int>&state, int maxNewSize = INT_MAX) {
		vector<vector<int>> out;
		for (int i = 0; i < (int)state.size(); ++i)
			for (int to = 0; to < state[i] && to <= maxNewSize; ++to) {
				vector<int> next = state; next[i] = to;
				out.push_back(canonicalPiles(next));
			}
		sort(out.begin(),out.end()); out.erase(unique(out.begin(),out.end()),out.end());
		return out;
	}

	// Moore Nim: reduce arbitrary positive amounts in 1..k distinct heaps.
	static vector<vector<int>> mooreMoves(const vector<int>&state, int k) {
		assert(k >= 1); vector<vector<int>> out; int n = state.size();
		function<void(int,int,vector<int>&)> dfs = [&](int at, int changed, vector<int>&next) {
			if (at == n) {
				if (changed) out.push_back(canonicalPiles(next));
				return;
			}
			dfs(at+1,changed,next);
			if (changed == k) return;
			int old = next[at];
			for (int x = 0; x < old; ++x) next[at] = x, dfs(at+1,changed+1,next);
			next[at] = old;
		};
		vector<int> next = state; dfs(0,0,next);
		sort(out.begin(),out.end()); out.erase(unique(out.begin(),out.end()),out.end());
		return out;
	}

	// Another Game: decrement any nonempty subset of positive heaps.
	static vector<vector<int>> anotherGameMoves(const vector<int>&state) {
		vector<vector<int>> out; int n = state.size();
		for (int mask = 1; mask < (1<<n); ++mask) {
			vector<int> to = state; bool legal = true;
			for (int i = 0; i < n; ++i) if (mask>>i&1) {
				if (!to[i]) { legal = false; break; }
				--to[i];
			}
			if (legal) out.push_back(canonicalPiles(to,false));
		}
		return out;
	}

	template<class State>
	static void printStates(const vector<State>&states, Solver<State>&solver,
			ostream&out = cout) {
		out << "state\tW/L\tSG\n";
		for (const State&state : states)
			out << state << '\t' << (solver.win(state) ? 'W' : 'L')
				<< '\t' << solver.grundy(state) << '\n';
	}

	static void printSequence(int maximum, Solver<int>&solver, ostream&out = cout) {
		out << "x\tW/L\tSG\tparity\n";
		for (int x = 0; x <= maximum; ++x)
			out << x << '\t' << (solver.win(x) ? 'W' : 'L') << '\t'
				<< solver.grundy(x) << '\t' << (x&1) << '\n';
	}

	static void printLosingPairs(int maximum, Solver<pair<int,int>>&solver,
			ostream&out = cout) {
		out << "losing pairs (a,b), a<=b\n";
		for (int a = 0; a <= maximum; ++a)
			for (int b = a; b <= maximum; ++b)
				if (!solver.win({a,b})) out << '(' << a << ',' << b << ") ";
		out << '\n';
	}

	static void printLosingTriples(int maximum, Solver<array<int,3>>&solver,
			ostream&out = cout) {
		out << "losing triples (a,b,c), a<=b<=c\n";
		for (int a = 0; a <= maximum; ++a)
			for (int b = a; b <= maximum; ++b)
				for (int c = b; c <= maximum; ++c)
					if (!solver.win({a,b,c})) out << '(' << a << ',' << b << ',' << c << ") ";
		out << '\n';
	}

	// Returns {preperiod,period}; {-1,-1} if the available suffix is insufficient.
	template<class T>
	static pair<int,int> findEventualPeriod(const vector<T>&a, int minimumRepeats = 3) {
		assert(minimumRepeats >= 2); int n = a.size();
		for (int pre = 0; pre < n; ++pre)
			for (int period = 1; pre + minimumRepeats*period <= n; ++period) {
				bool good = true;
				for (int i = pre+period; i < n; ++i)
					if (a[i] != a[i-period]) { good = false; break; }
				if (good) return {pre,period};
			}
		return {-1,-1};
	}

	// Return the first x where brute(x) != formula(x), or -1.
	template<class Brute, class Formula>
	static int firstCounterexample1D(int maximum, Brute brute, Formula formula) {
		for (int x = 0; x <= maximum; ++x) if (bool(brute(x)) != bool(formula(x))) return x;
		return -1;
	}

	// Return the first canonical pair, or {-1,-1}.
	template<class Brute, class Formula>
	static pair<int,int> firstCounterexample2D(int maximum, Brute brute, Formula formula) {
		for (int a = 0; a <= maximum; ++a)
			for (int b = a; b <= maximum; ++b)
				if (bool(brute(a,b)) != bool(formula(a,b))) return {a,b};
		return {-1,-1};
	}

	// Return the first canonical triple, or {-1,-1,-1}.
	template<class Brute, class Formula>
	static array<int,3> firstCounterexample3D(int maximum, Brute brute, Formula formula) {
		for (int a = 0; a <= maximum; ++a)
			for (int b = a; b <= maximum; ++b)
				for (int c = b; c <= maximum; ++c)
					if (bool(brute(a,b,c)) != bool(formula(a,b,c))) return {a,b,c};
		return {-1,-1,-1};
	}

	static void printIntegerAnnotations(int maximum, ostream&out = cout) {
		auto spf = GameTheory::smallestPrimeFactor(max(1,maximum));
		out << "x\tbin\tparity\tprime exponents\n";
		for (int x = 1; x <= maximum; ++x) {
			string bits;
			for (int y = x; y; y >>= 1) bits += char('0'+(y&1));
			reverse(bits.begin(),bits.end());
			out << x << '\t' << bits << '\t' << (x&1) << '\t';
			for (auto [p,e] : GameTheory::factorExponents(x,spf)) out << p << '^' << e << ' ';
			out << '\n';
		}
	}

	static void printPairAnnotations(int maximum, ostream&out = cout) {
		out << "a\tb\tdifference\txor\tsum\n";
		for (int a = 0; a <= maximum; ++a)
			for (int b = a; b <= maximum; ++b)
				out << a << '\t' << b << '\t' << b-a << '\t'
					<< (a^b) << '\t' << a+b << '\n';
	}

	// Exhaustively compare a formula with minimax over sorted pile tuples.
	template<class MoveGenerator, class Formula>
	static optional<vector<int>> checkPileFormula(int heaps, int maximum,
			MoveGenerator moves, Formula formula) {
		Solver<vector<int>> solver(moves); vector<int> state(heaps);
		optional<vector<int>> bad;
		function<void(int,int)> dfs = [&](int at, int low) {
			if (bad) return;
			if (at == heaps) {
				if (solver.win(state) != bool(formula(state))) bad = state;
				return;
			}
			for (int x = low; x <= maximum; ++x) state[at] = x, dfs(at+1,x);
		};
		dfs(0,0); return bad;
	}

	// Ready-made regression for the most commonly misquoted Moore-Nim rule.
	static optional<vector<int>> checkMooreFormula(int heaps, int maximum, int k) {
		return checkPileFormula(heaps,maximum,
			[k](const vector<int>&s) { return mooreMoves(s,k); },
			[k](const vector<int>&s) {
				vector<ll> piles(s.begin(),s.end());
				return GameTheory::mooreNimAtMostKWin(piles,k);
			});
	}
};

int main() {
	ios::sync_with_stdio(false);
	cin.tie(nullptr);

	// Example 1: inspect W/L and SG patterns for moves {1,3,4}.
	auto moves = GameTheoryLab::subtractionMoves({1,3,4});
	GameTheoryLab::Solver<int> game(moves);
	GameTheoryLab::printSequence(30,game);

	vector<int> win(31), sg(31);
	for (int x = 0; x <= 30; ++x) win[x] = game.win(x), sg[x] = game.grundy(x);
	auto [winPre,winPeriod] = GameTheoryLab::findEventualPeriod(win);
	auto [sgPre,sgPeriod] = GameTheoryLab::findEventualPeriod(sg);
	cout << "W/L candidate: preperiod=" << winPre << ", period=" << winPeriod << '\n';
	cout << "SG candidate:  preperiod=" << sgPre << ", period=" << sgPeriod << '\n';

	// Example 2: try to falsify the Moore-Nim formula on small states.
	auto bad = GameTheoryLab::checkMooreFormula(3,4,2);
	if (!bad) cout << "Moore formula passed all tested states\n";
	else {
		cout << "counterexample:";
		for (int pile : *bad) cout << ' ' << pile;
		cout << '\n';
	}
}
