# An O(N sqrt(N)) Dynamic Programming Optimization for Partition Numbers

This section presents an (O(N\sqrt{N})) optimization for a two-dimensional dynamic programming formulation of the integer partition problem. The method applies to dynamic programs indexed by two parameters, here denoted by `n` and `k`, where the transition structure follows a small number of regular geometric directions in the DP table. The key idea is to compute the low-`k` region explicitly, while aggregating the high-`k` region through line sums. This avoids the quadratic cost of evaluating every state independently.

Let `DP[n][k]` denote the number of partitions of `n` into exactly `k` positive parts. The classical recurrence for this quantity is:

```math
DP[n][k] = DP[n - 1][k - 1] + DP[n - k][k]
```

This recurrence has a simple combinatorial interpretation. The term `DP[n - 1][k - 1]` accounts for partitions in which one unit is removed in a way that decreases both the target sum and the number of parts. The term `DP[n - k][k]` accounts for partitions in which the number of parts remains equal to `k`, while one unit is removed from each of the `k` parts. These two cases cover all partitions of `n` into exactly `k` positive parts.

A direct algorithm evaluates every state `DP[n][k]` for the ranges:

```math
1 \le n \le N
```

```math
1 \le k \le N
```

This direct computation requires:

```math
O(N^2)
```

time, since the table has quadratic size. The optimization reduces this cost by separating the DP table into a small-`k` region and a large-`k` region. Define the threshold:

```math
B = \lfloor \sqrt{N} \rfloor
```

The small-`k` region consists of all states satisfying:

```math
k < B
```

This part contains `B` columns and `N` possible values of `n`. Therefore, it contains:

```math
O(NB)
```

states. Since `B = floor(sqrt(N))`, this cost becomes:

```math
O(N\sqrt{N})
```

The small-`k` region is computed directly using the original recurrence. These values are kept because the large-`k` recurrence reaches the boundary column:

```math
DP[n][B - 1]
```

This boundary appears when a transition from a large column `k` moves to column `k - 1`. In particular, when `k = B`, the transition reaches `B - 1`, which belongs to the already computed small-`k` region.

The large-`k` region consists of all states satisfying:

```math
k \ge B
```

Computing all states in this region one by one would still be too expensive. Instead, the method groups large-`k` states along structured lines in the DP table. Define the auxiliary value:

```math
x[s][t] = \sum_{k \ge B} DP[s - tk][k]
```

Here, `s` is a starting coordinate, and `t` is a slope parameter. For fixed `s` and `t`, the expression sums all large-`k` states located on a line of the form:

```math
DP[s - tk][k]
```

As `k` varies, the first coordinate changes linearly. This creates a line in the DP table. The value `x[s][t]` therefore represents a grouped line sum over the large-`k` region.

The large-`k` contribution to the partition number is obtained by setting `s = N` and `t = 0`:

```math
x[N][0] = \sum_{k \ge B} DP[N][k]
```

Thus, the full partition number `p(N)` is written as the sum of the direct small-`k` contribution and the grouped large-`k` contribution:

```math
p(N)
=
\sum_{1 \le k < B} DP[N][k]
+
x[N][0]
```

It remains to derive a recurrence for `x[s][t]`. Start from its definition:

```math
x[s][t]
=
\sum_{k \ge B} DP[s - tk][k]
```

Apply the original recurrence to each term in this sum:

```math
DP[n][k] = DP[n - 1][k - 1] + DP[n - k][k]
```

Substituting `n = s - tk` gives:

```math
x[s][t]
=
\sum_{k \ge B} DP[s - tk - 1][k - 1]
+
\sum_{k \ge B} DP[s - tk - k][k]
```

The second summation is transformed directly. Since:

```math
s - tk - k = s - (t + 1)k
```

we obtain:

```math
\sum_{k \ge B} DP[s - tk - k][k]
=
\sum_{k \ge B} DP[s - (t + 1)k][k]
=
x[s][t + 1]
```

This term corresponds to increasing the slope from `t` to `t + 1`.

The first summation requires a change of variable. Let:

```math
j = k - 1
```

Then:

```math
k = j + 1
```

The summand becomes:

```math
DP[s - t(j + 1) - 1][j]
```

which is equivalent to:

```math
DP[s - t - 1 - tj][j]
```

For all terms with:

```math
j \ge B
```

the summation matches the grouped value:

```math
x[s - t - 1][t]
```

One boundary term remains. This occurs when:

```math
j = B - 1
```

This term lies outside the large-`k` region, so it belongs to the precomputed small-`k` region. Its contribution is:

```math
DP[s - tB - 1][B - 1]
```

Combining the transformed first summation and second summation gives the recurrence:

```math
x[s][t]
=
x[s - t - 1][t]
+
DP[s - tB - 1][B - 1]
+
x[s][t + 1]
```

All invalid indices contribute zero. For example, if either condition holds:

```math
s - t - 1 < 0
```

```math
s - tB - 1 < 0
```

then the corresponding term is treated as zero.

This recurrence computes an entire grouped line sum rather than one DP state. It therefore compresses the computation of the large-`k` region. The evaluation order follows directly from the dependencies. For each fixed slope `t`, the values of `s` are processed from small to large. This ensures that the term:

```math
x[s - t - 1][t]
```

has already been computed, because it depends on a smaller value of `s`. The slopes are processed in decreasing order, which ensures that:

```math
x[s][t + 1]
```

is already available when computing `x[s][t]`.

The number of slope values that need to be processed is also bounded. A term in `x[s][t]` contributes only when:

```math
s - tk \ge 0
```

Since the large-`k` region satisfies:

```math
k \ge B
```

we have:

```math
tk \ge tB
```

For all relevant values `s <= N`, once `t` grows beyond approximately `sqrt(N)`, the expression `s - tk` becomes negative. Those states contribute zero. Therefore, only slopes up to the threshold are needed:

```math
0 \le t \le B
```

This gives:

```math
O(\sqrt{N})
```

possible slopes.

The complexity now follows from the table split. The small-`k` region has `B` columns and `N` rows, so it costs:

```math
O(NB)
```

The grouped large-`k` region processes `O(B)` slopes, and each slope processes `O(N)` values of `s`. Its cost is also:

```math
O(NB)
```

Since:

```math
B = \lfloor \sqrt{N} \rfloor
```

the total time complexity is:

```math
O(N\sqrt{N})
```

The memory usage stays linear. The small-`k` computation only needs the required DP columns and the boundary values:

```math
DP[n][B - 1]
```

The grouped recurrence only needs two consecutive slope layers:

```math
x[\cdot][t]
```

and:

```math
x[\cdot][t + 1]
```

Thus, the memory usage is:

```math
O(N)
```

The method replaces individual evaluation of high-`k` DP states with aggregated line sums. The original recurrence induces a recurrence over these line sums. This structure reduces the computation from:

```math
O(N^2)
```

to:

```math
O(N\sqrt{N})
```

while preserving linear memory.
