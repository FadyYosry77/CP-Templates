
/*
FUNCTION CATALOG
Function                         Use
---------------------------------------------------------------------------
getbit/setbit/clearbit           inspect or mutate one bit
togglebit/assignbit              flip or assign one bit
ispowerof2/atmostonebit          exact/optional single-bit checks
lowbit/removelowbit              isolate/remove lowest set bit
getleastbit/gethighestbit        extreme set-bit indexes
bitwidth/bitfloor/bitceil        representation and power-of-two boundaries
popcount/zerocount               count ones/zeros
firstbits/rangemask              construct low/interval masks
extractbits/replacebits          read/write a bit interval
containsmask/disjointmask        containment/intersection tests
issubmask/complementmask         subset relation and fixed-width complement
calc_ones/calc_zeros             bit frequency over [0,n]
calc_ones_range                  bit frequency over [l,r]
total_set_bits                   total popcount over an integer interval
bit_frequency                    per-bit counts for ranges or arrays
getsubmasks/getsupermasks        enumerate subset/superset masks
masks_with_k_bits               fixed-popcount masks via Gosper's hack
binary_to_gray/gray_to_binary    Gray-code conversion
reversebits/rotateleft/right     fixed-width permutation
hammingdistance/bitparity        comparison and parity
prefixxor/rangexor               period-4 XOR formulas
rangeand/rangeor                 bitwise operation over an integer interval
sum_pair_xor/and/or              unordered-pair contribution sums
sum_subarray_xor/and/or          contiguous-subarray contribution sums
count_distinct_subarray_or/and   compressed monotone bitwise states
count_subarray_xor               exact-XOR subarray count
unique_among_k                   bit-frequency modulo repetition count
two_unique_among_pairs           XOR lowbit partition
sum_subset_xor_mod               XOR sum over every subset
*/


/*
BIT MANIPULATION HELPERS -- C++17

All public helpers are independent free functions. Bit positions are 0-based.
Unless stated otherwise, masks and values are nonnegative signed long long,
so usable bit positions are 0..62.

QUICK USE
    Basic:
        getbit, setbit, clearbit, togglebit, assignbit
    Important:
        ispowerof2, lowbit, removelowbit
        getleastbit, gethighestbit, bitwidth, bitfloor, bitceil
    Counting:
        popcount, zerocount, calc_ones, calc_zeros
        calc_ones_range, total_set_bits, bit_frequency
    Masks:
        firstbits, rangemask, extractbits, replacebits
        containsmask, disjointmask, issubmask, complementmask
    Enumeration:
        setbitpositions, getsubmasks, getsupermasks
        next_combination_mask, masks_with_k_bits
    Transformations:
        binary_to_gray, gray_to_binary, reversebits
        rotateleft, rotateright, hammingdistance, bitparity
    Integer ranges:
        prefixxor, rangexor, rangeand, rangeor
    Contributions:
        sum_pair_xor, sum_pair_and, sum_pair_or
        sum_subarray_xor, sum_subarray_and, sum_subarray_or
    Common XOR tricks:
        count_subarray_xor, unique_among_pairs
        unique_among_k, unique_among_triples, two_unique_among_pairs

RECOGNITION
    One bit or mask mutation
        -> shifts, get/set/clear/toggle, firstbits/rangemask.
    Enumerate every subset of one mask
        -> sub=(sub-1)&mask; O(2^popcount(mask)).
    Enumerate every (mask,submask) pair over n bits
        -> nested submask loop; O(3^n).
    Enumerate only k selected positions
        -> Gosper's hack / masks_with_k_bits.
    XOR of every integer in [l,r]
        -> period-4 prefixxor.
    AND/OR of every integer in [l,r]
        -> only the common binary prefix survives; lower changing bits differ.
    Count bit b among integers 0..n
        -> blocks of 2^(b+1): zero half, then one half.
    Sum XOR/AND/OR over pairs
        -> process each bit: differ / both one / at least one one.
    Sum XOR over subarrays
        -> unordered pairs of prefix XOR values.
    Sum AND over subarrays
        -> runs of ones for each bit.
    Sum OR over subarrays
        -> all subarrays minus runs of zeros for each bit.
    Maximum subset XOR or XOR-space rank
        -> use an XOR basis, not these helpers.
    Maximum XOR against stored values
        -> use a binary trie.
    Sum over all submasks for every mask
        -> use SOS DP.
    XOR/AND/OR convolution
        -> use the corresponding FWT.

IMPORTANT IDENTITIES
    x & (x-1)    removes the lowest set bit.
    x & -x       isolates the lowest set bit (for positive signed ll).
    x | (x-1)    changes bits below the lowest set bit to ones.
    x & (x+1)    clears the trailing ones.
    x | (x+1)    sets the lowest zero bit.
    ~x & (x+1)   isolates the lowest zero bit in an unsigned/fixed-width view.
    x ^ (x>>1)   converts binary to Gray code.
    x ^ (x-1)    makes bits 0..lowest-set-bit all ones.

OVERFLOW
    1LL<<63 is never used. Width is at most 63 and the highest value bit is 62.
    Exact aggregate sums assert that their mathematical result fits ll.
*/

#include <bits/stdc++.h>
using namespace std;

using ll = long long;
using i128 = __int128_t;

// ============================================================
// BASIC BIT ACCESS
// ============================================================

// Problem:
//     Test one bit of a nonnegative signed-ll value.
// Returns:
//     true iff bit is set.
// Complexity:
//     O(1).
// Notes:
//     Requires 0<=bit<63.
bool getbit(ll x,int bit){
    assert(x>=0&&0<=bit&&bit<63);
    return (x>>bit)&1LL;
}

// Problem:
//     Turn one bit on.
// Returns:
//     x with bit set to one.
// Complexity:
//     O(1).
// Notes:
//     Requires x>=0 and 0<=bit<63.
ll setbit(ll x,int bit){
    assert(x>=0&&0<=bit&&bit<63);
    return x|(1LL<<bit);
}

// Problem:
//     Turn one bit off.
// Returns:
//     x with bit cleared.
// Complexity:
//     O(1).
// Notes:
//     Requires x>=0 and 0<=bit<63.
ll clearbit(ll x,int bit){
    assert(x>=0&&0<=bit&&bit<63);
    return x&~(1LL<<bit);
}

// Problem:
//     Flip one bit.
// Returns:
//     x with bit toggled.
// Complexity:
//     O(1).
// Notes:
//     Requires x>=0 and 0<=bit<63.
ll togglebit(ll x,int bit){
    assert(x>=0&&0<=bit&&bit<63);
    return x^(1LL<<bit);
}

// Problem:
//     Set one bit to a supplied Boolean value.
// Returns:
//     x with bit equal to value.
// Complexity:
//     O(1).
// Notes:
//     Reuses setbit/clearbit.
ll assignbit(ll x,int bit,bool value){
    return value?setbit(x,bit):clearbit(x,bit);
}

// ============================================================
// LOWEST / HIGHEST BIT AND POWERS OF TWO
// ============================================================

// Problem:
//     Check whether x is an exact positive power of two.
// Returns:
//     true iff x=2^k for some k>=0.
// Complexity:
//     O(1).
// Notes:
//     Zero and negative numbers return false.
bool ispowerof2(ll x){
    return x>0&&!(x&(x-1));
}

// Problem:
//     Find the least-significant set-bit position readably.
// Returns:
//     First set position, or -1 for x=0.
// Complexity:
//     O(63).
// Notes:
//     Requires x>=0.
int getleastbit(ll x){
    assert(x>=0);
    for(int bit=0;bit<63;bit++){
        if((x>>bit)&1LL)return bit;
    }
    return -1;
}

// Problem:
//     Find the least-significant set-bit position with a builtin.
// Returns:
//     First set position, or -1 for x=0.
// Complexity:
//     O(1).
// Notes:
//     Guards zero because __builtin_ctzll(0) is undefined.
int getleastbit_fast(ll x){
    assert(x>=0);
    return x?__builtin_ctzll((unsigned long long)x):-1;
}

// Problem:
//     Isolate the least-significant set bit.
// Returns:
//     x&-x; returns 0 for x=0.
// Complexity:
//     O(1).
// Notes:
//     Requires x>=0.
ll lowbit(ll x){
    assert(x>=0);
    return x&-x;
}

// Problem:
//     Remove the least-significant set bit.
// Returns:
//     x&(x-1); returns 0 for x=0.
// Complexity:
//     O(1).
// Notes:
//     Repeated use enumerates/removes all set bits.
ll removelowbit(ll x){
    assert(x>=0);
    return x?x&(x-1):0;
}

// Problem:
//     Find the most-significant set-bit position readably.
// Returns:
//     Highest set position, or -1 for x=0.
// Complexity:
//     O(63).
// Notes:
//     Requires x>=0.
int gethighestbit(ll x){
    assert(x>=0);
    for(int bit=62;bit>=0;bit--){
        if((x>>bit)&1LL)return bit;
    }
    return -1;
}

// Problem:
//     Find the most-significant set-bit position with a builtin.
// Returns:
//     Highest set position, or -1 for x=0.
// Complexity:
//     O(1).
// Notes:
//     Guards zero because __builtin_clzll(0) is undefined.
int gethighestbit_fast(ll x){
    assert(x>=0);
    return x?63-__builtin_clzll((unsigned long long)x):-1;
}

// Problem:
//     Count the minimum bits needed to represent x.
// Returns:
//     0 for x=0, otherwise highest-set-bit+1.
// Complexity:
//     O(1).
// Notes:
//     Requires x>=0.
int bitwidth(ll x){
    return gethighestbit_fast(x)+1;
}

// Problem:
//     Find the largest power of two not exceeding x.
// Returns:
//     0 for x=0, otherwise 2^floor(log2(x)).
// Complexity:
//     O(1).
// Notes:
//     C++17 replacement for std::bit_floor.
ll bitfloor(ll x){
    assert(x>=0);
    return x?1LL<<gethighestbit_fast(x):0;
}

// Problem:
//     Find the smallest power of two not smaller than x.
// Returns:
//     bitceil(0)=1; -1 if the answer exceeds 2^62.
// Complexity:
//     O(1).
// Notes:
//     C++17 replacement for std::bit_ceil on nonnegative signed ll.
ll bitceil(ll x){
    assert(x>=0);
    if(x<=1)return 1;
    int bit=gethighestbit_fast(x-1)+1;
    return bit<63?1LL<<bit:-1;
}

// Problem:
//     Find the smallest power of two strictly larger than x.
// Returns:
//     1 for x=0; -1 if the answer exceeds 2^62.
// Complexity:
//     O(1).
// Notes:
//     Unlike bitceil, a power-of-two input advances to the next power.
ll strictnextpowerof2(ll x){
    assert(x>=0);
    if(x==0)return 1;
    int bit=gethighestbit_fast(x)+1;
    return bit<63?1LL<<bit:-1;
}

// Problem:
//     Recover the exponent of an exact power of two.
// Returns:
//     k for x=2^k, otherwise -1.
// Complexity:
//     O(1).
// Notes:
//     Zero and negative inputs return -1.
int powerof2exponent(ll x){
    return ispowerof2(x)?getleastbit_fast(x):-1;
}

// Problem:
//     Check whether at most one bit is set.
// Returns:
//     true for zero and positive powers of two.
// Complexity:
//     O(1).
// Notes:
//     Differs from ispowerof2 because zero is accepted.
bool atmostonebit(ll x){
    return x>=0&&!(x&(x-1));
}

// Problem:
//     Check whether exactly k bits are set.
// Returns:
//     true iff popcount(x)=k.
// Complexity:
//     O(1).
// Notes:
//     Requires x>=0 and k>=0.
bool haskbits(ll x,int k){
    assert(x>=0&&k>=0);
    return __builtin_popcountll((unsigned long long)x)==k;
}

// ============================================================
// POPCOUNT AND FIXED-WIDTH BIT OPERATIONS
// ============================================================

// Problem:
//     Count set bits using a builtin.
// Returns:
//     Number of ones in the nonnegative signed-ll representation.
// Complexity:
//     O(1).
// Notes:
//     Requires x>=0.
int popcount(ll x){
    assert(x>=0);
    return __builtin_popcountll((unsigned long long)x);
}

// Problem:
//     Count set bits by repeatedly removing the lowest one.
// Returns:
//     Number of ones.
// Complexity:
//     O(number of set bits).
// Notes:
//     Educational alternative to popcount.
int popcount_manual(ll x){
    assert(x>=0);
    int answer=0;
    while(x){
        x&=x-1;
        answer++;
    }
    return answer;
}

// Problem:
//     Create a mask whose lowest k bits are ones.
// Returns:
//     0 for k=0 and LLONG_MAX for k=63.
// Complexity:
//     O(1).
// Notes:
//     Requires 0<=k<=63; avoids 1LL<<63.
ll firstbits(int k){
    assert(0<=k&&k<=63);
    if(k==63)return LLONG_MAX;
    return k?(1LL<<k)-1:0;
}

// Problem:
//     Count zero bits inside a selected low-bit width.
// Returns:
//     bits-popcount(x masked to the first bits positions).
// Complexity:
//     O(1).
// Notes:
//     Requires x>=0 and 0<=bits<=63.
int zerocount(ll x,int bits){
    assert(x>=0&&0<=bits&&bits<=63);
    return bits-popcount(x&firstbits(bits));
}

// Problem:
//     Create a mask with inclusive positions l..r set.
// Returns:
//     The requested fixed-width interval mask.
// Complexity:
//     O(1).
// Notes:
//     Requires 0<=l<=r<63.
ll rangemask(int l,int r){
    assert(0<=l&&l<=r&&r<63);
    return firstbits(r+1)^firstbits(l);
}

// Problem:
//     Extract an inclusive bit interval and move it to bit zero.
// Returns:
//     Bits x[l..r] shifted right by l.
// Complexity:
//     O(1).
// Notes:
//     Requires x>=0 and 0<=l<=r<63.
ll extractbits(ll x,int l,int r){
    assert(x>=0);
    return (x&rangemask(l,r))>>l;
}

// Problem:
//     Replace an inclusive bit interval.
// Returns:
//     x with [l,r] replaced by the low r-l+1 bits of value.
// Complexity:
//     O(1).
// Notes:
//     Higher bits of value are ignored.
ll replacebits(ll x,int l,int r,ll value){
    assert(x>=0&&value>=0);
    ll mask=rangemask(l,r);
    return (x&~mask)|((value&firstbits(r-l+1))<<l);
}

// Problem:
//     Check whether every bit of mask exists in x.
// Returns:
//     true iff (x&mask)==mask.
// Complexity:
//     O(1).
// Notes:
//     Requires nonnegative masks.
bool containsmask(ll x,ll mask){
    assert(x>=0&&mask>=0);
    return (x&mask)==mask;
}

// Problem:
//     Check whether two masks share no set bit.
// Returns:
//     true iff (x&mask)==0.
// Complexity:
//     O(1).
// Notes:
//     Requires nonnegative masks.
bool disjointmask(ll x,ll mask){
    assert(x>=0&&mask>=0);
    return !(x&mask);
}

// Problem:
//     Check a submask relationship.
// Returns:
//     true iff every bit of sub is also set in mask.
// Complexity:
//     O(1).
// Notes:
//     Equality and zero are valid submasks.
bool issubmask(ll sub,ll mask){
    assert(sub>=0&&mask>=0);
    return (sub&mask)==sub;
}

// Problem:
//     Complement a mask only inside its lowest bits positions.
// Returns:
//     (~mask)&firstbits(bits).
// Complexity:
//     O(1).
// Notes:
//     Requires mask to fit inside the selected width.
ll complementmask(ll mask,int bits){
    assert(mask>=0&&0<=bits&&bits<=63&&issubmask(mask,firstbits(bits)));
    return firstbits(bits)^mask;
}

// Problem:
//     Convert selected bit indexes into one mask.
// Returns:
//     OR of 1LL<<bit for all supplied positions.
// Complexity:
//     O(number of positions).
// Notes:
//     Duplicate positions are harmless; every bit must be in 0..62.
ll tomask(const vector<int>&bits){
    ll mask=0;
    for(int bit:bits)mask=setbit(mask,bit);
    return mask;
}

// Problem:
//     Return set-bit positions in increasing order.
// Returns:
//     Vector of every bit index whose bit is one.
// Complexity:
//     O(63).
// Notes:
//     Readable fixed-width implementation.
vector<int> setbitpositions(ll mask){
    assert(mask>=0);
    vector<int> answer;
    for(int bit=0;bit<63;bit++){
        if((mask>>bit)&1LL)answer.push_back(bit);
    }
    return answer;
}

// Problem:
//     Return set-bit positions using ctz.
// Returns:
//     Vector of every bit index in increasing order.
// Complexity:
//     O(number of set bits).
// Notes:
//     removelowbit guarantees ctz is never called on zero.
vector<int> setbitpositions_fast(ll mask){
    assert(mask>=0);
    vector<int> answer;
    while(mask){
        answer.push_back(getleastbit_fast(mask));
        mask=removelowbit(mask);
    }
    return answer;
}

// Problem:
//     Find the lowest zero position in a fixed width.
// Returns:
//     First zero bit, or -1 if all selected bits are one.
// Complexity:
//     O(bits).
// Notes:
//     Requires 0<=bits<=63.
int getleastzerobit(ll x,int bits=63){
    assert(x>=0&&0<=bits&&bits<=63);
    for(int bit=0;bit<bits;bit++){
        if(!((x>>bit)&1LL))return bit;
    }
    return -1;
}

// Problem:
//     Find the highest zero position in a fixed width.
// Returns:
//     Last zero bit, or -1 if all selected bits are one.
// Complexity:
//     O(bits).
// Notes:
//     Requires 0<=bits<=63.
int gethighestzerobit(ll x,int bits){
    assert(x>=0&&0<=bits&&bits<=63);
    for(int bit=bits-1;bit>=0;bit--){
        if(!((x>>bit)&1LL))return bit;
    }
    return -1;
}

// Problem:
//     Check whether x has form 2^k-1.
// Returns:
//     true for positive all-ones masks; false for zero.
// Complexity:
//     O(1).
// Notes:
//     Examples are 1,3,7,15.
bool isallones(ll x){
    return x>0&&(x==LLONG_MAX||!(x&(x+1)));
}

// Problem:
//     Check whether adjacent representation bits alternate.
// Returns:
//     true for positive values such as 1,2,5,10.
// Complexity:
//     O(1).
// Notes:
//     Leading zeros are ignored; zero returns false.
bool hasalternatingbits(ll x){
    if(x<=0)return false;
    ll adjacent=x^(x>>1);
    return isallones(adjacent);
}

// ============================================================
// BIT COUNTING OVER [0,N]
// ============================================================

// Problem:
//     Count integers x in [0,n] whose bit is set.
// Returns:
//     Number of ones at the 0-based bit position.
// Complexity:
//     O(1).
// Notes:
//     Returns 0 for n<0; uses i128 so bit=62 is safe.
ll calc_ones(ll n,int bit){
    assert(0<=bit&&bit<63);
    if(n<0)return 0;
    i128 count=(i128)n+1;
    i128 half=(i128)1<<bit;
    i128 cycle=half*2;
    i128 answer=count/cycle*half+max((i128)0,count%cycle-half);
    return (ll)answer;
}

// Problem:
//     Count integers x in [0,n] whose bit is zero.
// Returns:
//     (n+1)-calc_ones(n,bit), or 0 for n<0.
// Complexity:
//     O(1).
// Notes:
//     The interval includes zero.
ll calc_zeros(ll n,int bit){
    assert(0<=bit&&bit<63);
    if(n<0)return 0;
    return (ll)((i128)n+1-calc_ones(n,bit));
}

// Problem:
//     Count set occurrences at one bit in [l,r].
// Returns:
//     calc_ones(r,bit)-calc_ones(l-1,bit).
// Complexity:
//     O(1).
// Notes:
//     Requires 0<=l; reversed ranges return zero.
ll calc_ones_range(ll l,ll r,int bit){
    assert(l>=0);
    return l>r?0:calc_ones(r,bit)-calc_ones(l-1,bit);
}

// Problem:
//     Sum popcount(x) over every x in [0,n].
// Returns:
//     Exact total, or 0 for n<0.
// Complexity:
//     O(63).
// Notes:
//     Asserts if the mathematical sum does not fit ll.
ll total_set_bits(ll n){
    if(n<0)return 0;
    i128 answer=0;
    for(int bit=0;bit<63;bit++)answer+=calc_ones(n,bit);
    assert(answer<=LLONG_MAX);
    return (ll)answer;
}

// Problem:
//     Sum popcount(x) over every x in [l,r].
// Returns:
//     Prefix difference of total_set_bits.
// Complexity:
//     O(63).
// Notes:
//     Requires 0<=l and an exact result fitting ll.
ll total_set_bits(ll l,ll r){
    assert(l>=0);
    return l>r?0:total_set_bits(r)-total_set_bits(l-1);
}

// Problem:
//     Count every bit independently among integers [0,n].
// Returns:
//     answer[bit]=calc_ones(n,bit).
// Complexity:
//     O(bits).
// Notes:
//     Requires 0<=bits<=63.
vector<ll> bit_frequency(ll n,int bits=63){
    assert(0<=bits&&bits<=63);
    vector<ll> answer(bits);
    for(int bit=0;bit<bits;bit++)answer[bit]=calc_ones(n,bit);
    return answer;
}

// Problem:
//     Count every bit independently in an array.
// Returns:
//     answer[bit]=number of elements having that bit set.
// Complexity:
//     O(array size * bits).
// Notes:
//     Array values must be nonnegative and 0<=bits<=63.
vector<ll> bit_frequency(const vector<ll>&a,int bits=63){
    assert(0<=bits&&bits<=63);
    vector<ll> answer(bits);
    for(ll x:a){
        assert(x>=0);
        for(int bit=0;bit<bits;bit++)answer[bit]+=(x>>bit)&1LL;
    }
    return answer;
}

// ============================================================
// SUBMASK / SUPERMASK ENUMERATION
// ============================================================

// Problem:
//     Enumerate every submask, including mask and zero.
// Returns:
//     Submasks in decreasing numeric order.
// Complexity:
//     O(2^popcount(mask)) time and output space.
// Notes:
//     The break after zero avoids unsigned-style wraparound loops.
vector<ll> getsubmasks(ll mask){
    assert(mask>=0);
    vector<ll> answer;
    for(ll sub=mask;;sub=(sub-1)&mask){
        answer.push_back(sub);
        if(sub==0)break;
    }
    return answer;
}

// Problem:
//     Enumerate proper nonzero submasks.
// Returns:
//     All submasks except mask itself and zero.
// Complexity:
//     O(2^popcount(mask)).
// Notes:
//     Returns empty for masks with at most one set bit.
vector<ll> getpropernonzerosubmasks(ll mask){
    assert(mask>=0);
    vector<ll> answer;
    if(!mask)return answer;
    for(ll sub=(mask-1)&mask;sub;sub=(sub-1)&mask)answer.push_back(sub);
    return answer;
}

// Problem:
//     Enumerate supermasks of mask inside an n-bit universe.
// Returns:
//     Every sup with mask subset of sup<2^n.
// Complexity:
//     O(2^(n-popcount(mask))) time and output space.
// Notes:
//     Requires 0<=n<=63 and mask to fit n bits.
vector<ll> getsupermasks(ll mask,int n){
    assert(0<=n&&n<=63&&mask>=0&&issubmask(mask,firstbits(n)));
    ll freebits=complementmask(mask,n);
    vector<ll> answer;
    for(ll sub=freebits;;sub=(sub-1)&freebits){
        answer.push_back(mask|sub);
        if(sub==0)break;
    }
    return answer;
}

// Problem:
//     Find the next larger mask having the same popcount.
// Returns:
//     Gosper successor, or -1 when none fits nonnegative ll.
// Complexity:
//     O(1).
// Notes:
//     Requires x>0 for a successor candidate.
ll next_combination_mask(ll x){
    if(x<=0)return -1;
    unsigned long long value=(unsigned long long)x;
    unsigned long long smallest=value&(~value+1);
    unsigned long long ripple=value+smallest;
    if(ripple>(unsigned long long)LLONG_MAX)return -1;
    unsigned long long ones=((value^ripple)/smallest)>>2;
    return (ll)(ripple|ones);
}

// Problem:
//     Enumerate all n-bit masks with exactly k set bits.
// Returns:
//     Masks in increasing numeric order.
// Complexity:
//     O(C(n,k)) time and output space.
// Notes:
//     Uses Gosper's hack; supports 0<=n<=63.
vector<ll> masks_with_k_bits(int n,int k){
    assert(0<=n&&n<=63&&0<=k&&k<=n);
    if(k==0)return {0};
    unsigned long long limit=1ULL<<n;
    unsigned long long mask=(1ULL<<k)-1;
    vector<ll> answer;
    while(mask<limit){
        answer.push_back((ll)mask);
        ll next=next_combination_mask((ll)mask);
        if(next<0)break;
        mask=(unsigned long long)next;
    }
    return answer;
}

// ============================================================
// GRAY CODE, REVERSAL, ROTATION AND DISTANCE
// ============================================================

// Problem:
//     Convert a nonnegative binary integer to Gray code.
// Returns:
//     x^(x>>1).
// Complexity:
//     O(1).
// Notes:
//     Consecutive binary indexes have Gray codes differing in one bit.
ll binary_to_gray(ll x){
    assert(x>=0);
    return x^(x>>1);
}

// Problem:
//     Convert Gray code back to binary.
// Returns:
//     The original nonnegative binary integer.
// Complexity:
//     O(log value).
// Notes:
//     Repeatedly XORs all shifted Gray prefixes.
ll gray_to_binary(ll gray){
    assert(gray>=0);
    ll answer=0;
    for(;gray;gray>>=1)answer^=gray;
    return answer;
}

// Problem:
//     Compute Hamming distance between two bit patterns.
// Returns:
//     popcount(a^b).
// Complexity:
//     O(1).
// Notes:
//     Requires nonnegative values.
int hammingdistance(ll a,ll b){
    assert(a>=0&&b>=0);
    return popcount(a^b);
}

// Problem:
//     Compute parity of the set-bit count.
// Returns:
//     0 for even popcount, 1 for odd.
// Complexity:
//     O(1).
// Notes:
//     Uses the guarded nonnegative signed-ll domain.
int bitparity(ll x){
    assert(x>=0);
    return __builtin_parityll((unsigned long long)x);
}

// Problem:
//     Reverse the lowest bits positions.
// Returns:
//     Fixed-width bit reversal.
// Complexity:
//     O(bits).
// Notes:
//     Bits above the selected width are ignored.
ll reversebits(ll x,int bits){
    assert(x>=0&&0<=bits&&bits<=63);
    ll answer=0;
    for(int bit=0;bit<bits;bit++){
        if((x>>bit)&1LL)answer|=1LL<<(bits-1-bit);
    }
    return answer;
}

// Problem:
//     Rotate a fixed-width value left.
// Returns:
//     Lowest bits positions rotated by shift.
// Complexity:
//     O(1).
// Notes:
//     Requires 1<=bits<=63; negative shifts are normalized.
ll rotateleft(ll x,int shift,int bits){
    assert(x>=0&&1<=bits&&bits<=63);
    unsigned long long mask=(unsigned long long)firstbits(bits);
    unsigned long long value=(unsigned long long)x&mask;
    shift%=bits;
    if(shift<0)shift+=bits;
    if(!shift)return (ll)value;
    return (ll)(((value<<shift)&mask)|(value>>(bits-shift)));
}

// Problem:
//     Rotate a fixed-width value right.
// Returns:
//     Lowest bits positions rotated by shift.
// Complexity:
//     O(1).
// Notes:
//     Implemented through the equivalent left rotation.
ll rotateright(ll x,int shift,int bits){
    assert(1<=bits&&bits<=63);
    shift%=bits;
    if(shift<0)shift+=bits;
    return rotateleft(x,bits-shift,bits);
}

// Problem:
//     Find the highest bit at which x and y differ.
// Returns:
//     gethighestbit(x^y), or -1 if equal.
// Complexity:
//     O(1).
// Notes:
//     Useful for binary lexicographic/common-prefix reasoning.
int highestdifferentbit(ll x,ll y){
    assert(x>=0&&y>=0);
    return gethighestbit_fast(x^y);
}

// Problem:
//     Find the lowest bit at which x and y differ.
// Returns:
//     getleastbit(x^y), or -1 if equal.
// Complexity:
//     O(1).
// Notes:
//     Useful for parity-class partitions.
int lowestdifferentbit(ll x,ll y){
    assert(x>=0&&y>=0);
    return getleastbit_fast(x^y);
}

// Problem:
//     Count equal leading bits inside a fixed width.
// Returns:
//     Length of the common high-bit prefix.
// Complexity:
//     O(1).
// Notes:
//     Both values must fit inside the chosen width.
int commonprefixbits(ll x,ll y,int bits){
    assert(x>=0&&y>=0&&0<=bits&&bits<=63);
    assert(issubmask(x,firstbits(bits))&&issubmask(y,firstbits(bits)));
    int different=highestdifferentbit(x,y);
    return different<0?bits:bits-different-1;
}

// ============================================================
// XOR / AND / OR OVER INTEGER RANGES
// ============================================================

// Problem:
//     XOR every integer in [0,n].
// Returns:
//     Period-4 prefix XOR; returns 0 for n<0.
// Complexity:
//     O(1).
// Notes:
//     XOR over [1,n] is numerically identical because zero changes nothing.
ll prefixxor(ll n){
    if(n<0)return 0;
    if(n%4==0)return n;
    if(n%4==1)return 1;
    if(n%4==2)return n+1;
    return 0;
}

// Problem:
//     XOR every integer in [l,r].
// Returns:
//     prefixxor(r)^prefixxor(l-1).
// Complexity:
//     O(1).
// Notes:
//     Requires 0<=l; reversed ranges return zero.
ll rangexor(ll l,ll r){
    assert(l>=0);
    return l>r?0:prefixxor(r)^prefixxor(l-1);
}

// Problem:
//     AND every integer in [l,r].
// Returns:
//     The common binary prefix followed by zeros.
// Complexity:
//     O(log(max(l,r))).
// Notes:
//     Requires 0<=l<=r.
ll rangeand(ll l,ll r){
    assert(0<=l&&l<=r);
    int shifts=0;
    while(l!=r){
        l>>=1;
        r>>=1;
        shifts++;
    }
    return shifts==63?0:l<<shifts;
}

// Problem:
//     OR every integer in [l,r].
// Returns:
//     Shared high prefix followed by ones at every changing lower bit.
// Complexity:
//     O(1).
// Notes:
//     Requires 0<=l<=r.
ll rangeor(ll l,ll r){
    assert(0<=l&&l<=r);
    if(l==r)return l;
    int different=gethighestbit_fast(l^r);
    return (l&~firstbits(different+1))|firstbits(different+1);
}

// ============================================================
// BIT CONTRIBUTION FUNCTIONS
// ============================================================

// Problem:
//     Sum a[i]^a[j] over unordered pairs i<j.
// Returns:
//     Exact XOR-pair contribution sum.
// Complexity:
//     O(63*n).
// Notes:
//     At each bit, ones*zeros pairs contribute; result must fit ll.
ll sum_pair_xor(const vector<ll>&a){
    vector<ll> ones=bit_frequency(a);
    i128 answer=0;
    for(int bit=0;bit<63;bit++){
        answer+=(i128)ones[bit]*(a.size()-ones[bit])*(1LL<<bit);
    }
    assert(answer<=LLONG_MAX);
    return (ll)answer;
}

// Problem:
//     Sum a[i]&a[j] over unordered pairs i<j.
// Returns:
//     Exact AND-pair contribution sum.
// Complexity:
//     O(63*n).
// Notes:
//     C(ones,2) pairs contribute at each bit; result must fit ll.
ll sum_pair_and(const vector<ll>&a){
    vector<ll> ones=bit_frequency(a);
    i128 answer=0;
    for(int bit=0;bit<63;bit++){
        answer+=(i128)ones[bit]*(ones[bit]-1)/2*(1LL<<bit);
    }
    assert(answer<=LLONG_MAX);
    return (ll)answer;
}

// Problem:
//     Sum a[i]|a[j] over unordered pairs i<j.
// Returns:
//     Exact OR-pair contribution sum.
// Complexity:
//     O(63*n).
// Notes:
//     Total pairs minus C(zeros,2) contribute at each bit.
ll sum_pair_or(const vector<ll>&a){
    vector<ll> ones=bit_frequency(a);
    i128 total=(i128)a.size()*(a.size()-1)/2;
    i128 answer=0;
    for(int bit=0;bit<63;bit++){
        i128 zeros=(ll)a.size()-ones[bit];
        answer+=(total-zeros*(zeros-1)/2)*(1LL<<bit);
    }
    assert(answer<=LLONG_MAX);
    return (ll)answer;
}

// Problem:
//     Sum pairwise XOR modulo an arbitrary positive modulus.
// Returns:
//     sum_pair_xor reduced modulo mod without exact-sum overflow.
// Complexity:
//     O(63*n).
// Notes:
//     Uses i128 for each modular multiplication.
ll sum_pair_xor_mod(const vector<ll>&a,ll mod){
    assert(mod>0);
    vector<ll> ones=bit_frequency(a);
    ll answer=0;
    ll power=1%mod;
    for(int bit=0;bit<63;bit++){
        answer=(ll)((answer+(i128)(ones[bit]%mod)*
                (((ll)a.size()-ones[bit])%mod)%mod*power)%mod);
        power=(ll)((i128)power*2%mod);
    }
    return answer;
}

// Problem:
//     Sum XOR over all nonempty contiguous subarrays.
// Returns:
//     Sum over XOR of every prefix-XOR pair.
// Complexity:
//     O(63*n).
// Notes:
//     Includes prefix 0; exact answer must fit ll.
ll sum_subarray_xor(const vector<ll>&a){
    vector<ll> ones(63);
    ll prefix=0;
    for(ll x:a){
        assert(x>=0);
        prefix^=x;
        for(int bit=0;bit<63;bit++)ones[bit]+=(prefix>>bit)&1LL;
    }
    i128 answer=0;
    ll prefixes=(ll)a.size()+1;
    for(int bit=0;bit<63;bit++){
        answer+=(i128)ones[bit]*(prefixes-ones[bit])*(1LL<<bit);
    }
    assert(answer<=LLONG_MAX);
    return (ll)answer;
}

// Problem:
//     Sum AND over all nonempty contiguous subarrays.
// Returns:
//     Per-bit sum of len*(len+1)/2 over runs of ones.
// Complexity:
//     O(63*n).
// Notes:
//     Exact answer must fit ll.
ll sum_subarray_and(const vector<ll>&a){
    i128 answer=0;
    for(int bit=0;bit<63;bit++){
        ll run=0;
        for(ll x:a){
            assert(x>=0);
            if((x>>bit)&1LL)run++;
            else answer+=(i128)run*(run+1)/2*(1LL<<bit),run=0;
        }
        answer+=(i128)run*(run+1)/2*(1LL<<bit);
    }
    assert(answer<=LLONG_MAX);
    return (ll)answer;
}

// Problem:
//     Sum OR over all nonempty contiguous subarrays.
// Returns:
//     For each bit, total subarrays minus all-zero runs.
// Complexity:
//     O(63*n).
// Notes:
//     Exact answer must fit ll.
ll sum_subarray_or(const vector<ll>&a){
    i128 total=(i128)a.size()*(a.size()+1)/2;
    i128 answer=0;
    for(int bit=0;bit<63;bit++){
        i128 zero_subarrays=0;
        ll run=0;
        for(ll x:a){
            assert(x>=0);
            if(!((x>>bit)&1LL))run++;
            else zero_subarrays+=(i128)run*(run+1)/2,run=0;
        }
        zero_subarrays+=(i128)run*(run+1)/2;
        answer+=(total-zero_subarrays)*(1LL<<bit);
    }
    assert(answer<=LLONG_MAX);
    return (ll)answer;
}

// Problem:
//     Count distinct bitwise-OR values of all nonempty subarrays.
// Returns:
//     Number of distinct results.
// Complexity:
//     O(n*63*log(answer size)).
// Notes:
//     OR states ending at one position compress to O(63).
ll count_distinct_subarray_or(const vector<ll>&a){
    set<ll> all;
    vector<ll> current;
    for(ll x:a){
        assert(x>=0);
        vector<ll> next={x};
        for(ll value:current){
            ll merged=value|x;
            if(next.back()!=merged)next.push_back(merged);
        }
        current.swap(next);
        for(ll value:current)all.insert(value);
    }
    return (ll)all.size();
}

// Problem:
//     Count distinct bitwise-AND values of all nonempty subarrays.
// Returns:
//     Number of distinct results.
// Complexity:
//     O(n*63*log(answer size)).
// Notes:
//     AND states ending at one position compress because bits only disappear.
ll count_distinct_subarray_and(const vector<ll>&a){
    set<ll> all;
    vector<ll> current;
    for(ll x:a){
        assert(x>=0);
        vector<ll> next={x};
        for(ll value:current){
            ll merged=value&x;
            if(next.back()!=merged)next.push_back(merged);
        }
        current.swap(next);
        for(ll value:current)all.insert(value);
    }
    return (ll)all.size();
}

// ============================================================
// COMMON XOR TRICKS
// ============================================================

// Problem:
//     Count subarrays whose XOR equals target.
// Returns:
//     Number of pairs of prefixes differing by target.
// Complexity:
//     Expected O(n).
// Notes:
//     Uses an unordered_map of previous prefix-XOR frequencies.
ll count_subarray_xor(const vector<ll>&a,ll target){
    assert(target>=0);
    unordered_map<ll,ll> frequency;
    frequency.reserve(a.size()*2+1);
    frequency[0]=1;
    ll prefix=0,answer=0;
    for(ll x:a){
        assert(x>=0);
        prefix^=x;
        answer+=frequency[prefix^target];
        frequency[prefix]++;
    }
    return answer;
}

// Problem:
//     Recover the unique value when every other value occurs an even count.
// Returns:
//     XOR of all elements.
// Complexity:
//     O(n).
// Notes:
//     The usual statement has every other value appearing exactly twice.
ll unique_among_pairs(const vector<ll>&a){
    ll answer=0;
    for(ll x:a)answer^=x;
    return answer;
}

// Problem:
//     Recover one value occurring once while every other value occurs k times.
// Returns:
//     The exceptional nonnegative value.
// Complexity:
//     O(63*n).
// Notes:
//     Requires k>=2 and exactly the stated frequency pattern.
ll unique_among_k(const vector<ll>&a,int k){
    assert(k>=2);
    ll answer=0;
    for(int bit=0;bit<63;bit++){
        ll count=0;
        for(ll x:a){
            assert(x>=0);
            count+=(x>>bit)&1LL;
        }
        if(count%k)answer|=1LL<<bit;
    }
    return answer;
}

// Problem:
//     Recover one value occurring once while every other value occurs three times.
// Returns:
//     The exceptional nonnegative value.
// Complexity:
//     O(63*n).
// Notes:
//     Thin readable specialization of unique_among_k.
ll unique_among_triples(const vector<ll>&a){
    return unique_among_k(a,3);
}

// Problem:
//     Recover two unique values while every other value occurs in pairs.
// Returns:
//     The two unique values in increasing order.
// Complexity:
//     O(n).
// Notes:
//     Their XOR supplies a distinguishing low bit; assumes exactly two uniques.
pair<ll,ll> two_unique_among_pairs(const vector<ll>&a){
    ll xr=0;
    for(ll x:a)xr^=x;
    assert(xr!=0);
    ll split=lowbit(xr);
    ll first=0,second=0;
    for(ll x:a){
        if(x&split)first^=x;
        else second^=x;
    }
    if(first>second)swap(first,second);
    return {first,second};
}

// Problem:
//     Compute 2^exponent modulo mod.
// Returns:
//     Modular power used by subset-XOR formulas.
// Complexity:
//     O(log exponent).
// Notes:
//     Requires exponent>=0 and mod>0.
ll power2mod(int exponent,ll mod){
    assert(exponent>=0&&mod>0);
    ll base=2%mod,answer=1%mod;
    while(exponent){
        if(exponent&1)answer=(ll)((i128)answer*base%mod);
        base=(ll)((i128)base*base%mod);
        exponent>>=1;
    }
    return answer;
}

// Problem:
//     Sum XOR-values over all subsets modulo mod.
// Returns:
//     OR(all values)*2^(n-1) modulo mod for nonempty arrays.
// Complexity:
//     O(n+log n).
// Notes:
//     Includes the empty subset, whose XOR is zero.
ll sum_subset_xor_mod(const vector<ll>&a,ll mod){
    assert(mod>0);
    if(a.empty())return 0;
    ll support=0;
    for(ll x:a){
        assert(x>=0);
        support|=x;
    }
    return (ll)((i128)(support%mod)*power2mod((int)a.size()-1,mod)%mod);
}

#ifdef LOCAL
namespace bit_tests {

// Problem:
//     Brute-count one bit over [0,n] for differential testing.
// Returns:
//     Directly enumerated count.
// Complexity:
//     O(n).
// Notes:
//     LOCAL-only oracle for calc_ones.
ll brute_calc_ones(ll n,int bit){
    ll answer=0;
    for(ll x=0;x<=n;x++)answer+=(x>>bit)&1LL;
    return answer;
}

// Problem:
//     Brute-force AND over an integer interval.
// Returns:
//     l&(l+1)&...&r.
// Complexity:
//     O(r-l+1).
// Notes:
//     LOCAL-only oracle for rangeand.
ll brute_rangeand(ll l,ll r){
    ll answer=LLONG_MAX;
    for(ll x=l;x<=r;x++)answer&=x;
    return answer;
}

// Problem:
//     Brute-force OR over an integer interval.
// Returns:
//     l|(l+1)|...|r.
// Complexity:
//     O(r-l+1).
// Notes:
//     LOCAL-only oracle for rangeor.
ll brute_rangeor(ll l,ll r){
    ll answer=0;
    for(ll x=l;x<=r;x++)answer|=x;
    return answer;
}

// Problem:
//     Validate scalar, power, Gray, rotation and reversal helpers.
// Returns:
//     Nothing; asserts on the first mismatch.
// Complexity:
//     Small exhaustive LOCAL test.
// Notes:
//     Includes zero and every signed-ll power of two.
void test_scalar_helpers(){
    for(ll x=0;x<=10000;x++){
        assert(getleastbit(x)==getleastbit_fast(x));
        assert(gethighestbit(x)==gethighestbit_fast(x));
        assert(popcount(x)==popcount_manual(x));
        assert(bitwidth(x)==(x?gethighestbit(x)+1:0));
        assert(bitfloor(x)==(x?1LL<<gethighestbit(x):0));
        assert(binary_to_gray(gray_to_binary(binary_to_gray(x)))==binary_to_gray(x));
        assert(gray_to_binary(binary_to_gray(x))==x);
        for(int bits=1;bits<=16;bits++){
            ll masked=x&firstbits(bits);
            assert(rotateright(rotateleft(masked,7,bits),7,bits)==masked);
            assert(reversebits(reversebits(masked,bits),bits)==masked);
        }
    }

    for(int exponent=0;exponent<=62;exponent++){
        ll value=1LL<<exponent;
        assert(ispowerof2(value));
        assert(powerof2exponent(value)==exponent);
        assert(bitceil(value)==value);
        assert(bitfloor(value)==value);
    }

    for(ll x=0;x<=100000;x++){
        ll ceil=1;
        while(ceil<x)ceil*=2;
        assert(bitceil(x)==ceil);
    }
}

// Problem:
//     Validate periodic bit counts and integer-range operations.
// Returns:
//     Nothing; asserts on the first mismatch.
// Complexity:
//     Small exhaustive LOCAL test.
// Notes:
//     Prints the first calc_ones counterexample before failing.
void test_counting(){
    for(ll n=0;n<=300;n++){
        ll total=0;
        for(ll x=0;x<=n;x++)total+=popcount(x);
        assert(total_set_bits(n)==total);
        for(int bit=0;bit<12;bit++){
            if(calc_ones(n,bit)!=brute_calc_ones(n,bit)){
                cerr<<"failed n = "<<n<<" bit = "<<bit<<'\n';
                assert(false);
            }
            assert(calc_zeros(n,bit)==n+1-brute_calc_ones(n,bit));
        }
    }

    for(ll l=0;l<=100;l++){
        for(ll r=l;r<=100;r++){
            for(int bit=0;bit<10;bit++){
                ll answer=0;
                for(ll x=l;x<=r;x++)answer+=(x>>bit)&1LL;
                assert(calc_ones_range(l,r,bit)==answer);
            }
            assert(rangeand(l,r)==brute_rangeand(l,r));
            assert(rangeor(l,r)==brute_rangeor(l,r));
            ll xr=0;
            for(ll x=l;x<=r;x++)xr^=x;
            assert(rangexor(l,r)==xr);
        }
    }
}

// Problem:
//     Validate submask, supermask and fixed-popcount enumeration.
// Returns:
//     Nothing; asserts on duplicates, omissions or invalid masks.
// Complexity:
//     O(3^10) for the configured LOCAL range.
// Notes:
//     Compares Gosper output with brute popcount filtering.
void test_masks(){
    for(int n=0;n<=10;n++){
        ll limit=1LL<<n;
        for(ll mask=0;mask<limit;mask++){
            vector<ll> sub=getsubmasks(mask);
            set<ll> unique(sub.begin(),sub.end());
            assert(unique.size()==sub.size());
            assert((ll)sub.size()==(1LL<<popcount(mask)));
            for(ll value:sub)assert(issubmask(value,mask));

            vector<ll> super=getsupermasks(mask,n);
            set<ll> unique_super(super.begin(),super.end());
            assert(unique_super.size()==super.size());
            assert((ll)super.size()==(1LL<<(n-popcount(mask))));
            for(ll value:super)assert(issubmask(mask,value));
        }

        for(int k=0;k<=n;k++){
            vector<ll> actual=masks_with_k_bits(n,k);
            vector<ll> expected;
            for(ll mask=0;mask<limit;mask++){
                if(popcount(mask)==k)expected.push_back(mask);
            }
            assert(actual==expected);
        }
    }
}

// Problem:
//     Validate pair/subarray contribution formulas against nested loops.
// Returns:
//     Nothing; asserts on the first mismatch.
// Complexity:
//     Randomized small-array LOCAL test.
// Notes:
//     Also validates compressed distinct subarray OR/AND states.
void test_contributions(){
    mt19937_64 rng(123456789);

    for(int test=0;test<5000;test++){
        int n=(int)(rng()%10);
        vector<ll> a(n);
        for(ll&x:a)x=rng()%128;

        ll pair_xor=0,pair_and=0,pair_or=0;
        for(int i=0;i<n;i++){
            for(int j=i+1;j<n;j++){
                pair_xor+=a[i]^a[j];
                pair_and+=a[i]&a[j];
                pair_or+=a[i]|a[j];
            }
        }

        assert(sum_pair_xor(a)==pair_xor);
        assert(sum_pair_and(a)==pair_and);
        assert(sum_pair_or(a)==pair_or);

        ll sub_xor=0,sub_and=0,sub_or=0;
        set<ll> distinct_or,distinct_and;
        unordered_map<ll,ll> xor_count;

        for(int l=0;l<n;l++){
            ll xr=0,an=LLONG_MAX,orr=0;
            for(int r=l;r<n;r++){
                xr^=a[r];
                an&=a[r];
                orr|=a[r];
                sub_xor+=xr;
                sub_and+=an;
                sub_or+=orr;
                distinct_or.insert(orr);
                distinct_and.insert(an);
                xor_count[xr]++;
            }
        }

        assert(sum_subarray_xor(a)==sub_xor);
        assert(sum_subarray_and(a)==sub_and);
        assert(sum_subarray_or(a)==sub_or);
        assert(count_distinct_subarray_or(a)==(ll)distinct_or.size());
        assert(count_distinct_subarray_and(a)==(ll)distinct_and.size());

        for(ll target=0;target<32;target++){
            assert(count_subarray_xor(a,target)==xor_count[target]);
        }
    }
}

// Problem:
//     Validate repeated-frequency XOR recovery tricks.
// Returns:
//     Nothing; asserts on the first mismatch.
// Complexity:
//     O(number of generated test values * 63).
// Notes:
//     Covers one-among-triples and two-among-pairs.
void test_unique_helpers(){
    for(int test=0;test<1000;test++){
        vector<ll> triples;
        ll unique=test%101;
        triples.push_back(unique);
        for(int value=101;value<106;value++){
            for(int repeat=0;repeat<3;repeat++)triples.push_back(value);
        }
        assert(unique_among_triples(triples)==unique);

        vector<ll> pairs={3,3,8,10,10,17};
        assert(two_unique_among_pairs(pairs)==make_pair(8LL,17LL));
    }
}

// Problem:
//     Run every LOCAL bit-helper test group.
// Returns:
//     Nothing; prints success after all assertions pass.
// Complexity:
//     Intended only for local validation.
// Notes:
//     Compile with -DLOCAL.
void run(){
    test_scalar_helpers();
    test_counting();
    test_masks();
    test_contributions();
    test_unique_helpers();
    cerr<<"All bit-helper tests passed.\n";
}

}

// Problem:
//     Execute the LOCAL stress suite.
// Returns:
//     Process exit status zero on success.
// Complexity:
//     Sum of the LOCAL test groups.
// Notes:
//     This main exists only when LOCAL is defined.
int main(){
    bit_tests::run();
}
#endif
