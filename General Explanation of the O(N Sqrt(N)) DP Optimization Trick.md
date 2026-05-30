# The O(N sqrt(N)) DP Trick for Partition Numbers

Use this optimization when your DP table has two coordinates, `n` and `k`, and the transition moves in two regular directions.

For the partition DP, define `DP[n][k]` as the number of partitions of `n` into exactly `k` positive parts.

The standard recurrence is:

```math
DP[n][k] = DP[n - 1][k - 1] + DP[n - k][k]
```

The first term, `DP[n - 1][k - 1]`, reduces each part by `1`.

The second term, `DP[n - k][k]`, keeps the number of parts equal to `k`, and removes `1` from each of the `k` parts.

This recurrence splits the partitions into two cases. In the first case, each part contributes at least `1`, so removing that contribution reduces both the sum and the number of parts. In the second case, the number of parts stays fixed, and one unit is removed from every part.

A direct solution fills all states `DP[n][k]` for the following range:

```math
1 \le n \le N
```

```math
1 \le k \le N
```

This costs:

```math
O(N^2)
```

The goal is to avoid computing all large-`k` states one by one.

The optimization starts by splitting the DP table at:

```math
B = \lfloor \sqrt{N} \rfloor
```

For small values of `k`, meaning:

```math
k < B
```

compute the DP normally.

This region has `B` columns and `N` rows, so its cost is:

```math
O(NB) = O(N\sqrt{N})
```

This small-`k` region gives all values `DP[n][k]` for:

```math
k < B
```

Keep these values because the large-`k` computation touches their boundary.

The important boundary column is:

```math
DP[n][B - 1]
```

This boundary appears when a transition from a large column `k` moves to column `k - 1`. When `k = B`, the transition reaches `B - 1`, which belongs to the small-`k` region.

For large values of `k`, meaning:

```math
k \ge B
```

avoid computing every state separately.

Group many states together along lines.

Define:

```math
x[s][t] = \sum_{k \ge B} DP[s - tk][k]
```

Here, `t` represents the slope of the line.

For fixed `s` and `t`, this expression sums all large-`k` DP values lying on one line in the DP table.

Each term in this line has the form:

```math
DP[s - tk][k]
```

As `k` changes, the first coordinate changes linearly. This structure allows many DP states to be handled together.

The value needed from the large-`k` region is:

```math
x[N][0] = \sum_{k \ge B} DP[N][k]
```

So the final answer becomes:

```math
p(N)
=
\sum_{1 \le k < B} DP[N][k]
+
x[N][0]
```

The first term is the contribution of the small-`k` region.

The second term is the contribution of the large-`k` region.

Now apply the original recurrence to every term inside `x[s][t]`.

Start with:

```math
x[s][t]
=
\sum_{k \ge B} DP[s - tk][k]
```

Using:

```math
DP[n][k] = DP[n - 1][k - 1] + DP[n - k][k]
```

we get:

```math
x[s][t]
=
\sum_{k \ge B} DP[s - tk - 1][k - 1]
+
\sum_{k \ge B} DP[s - tk - k][k]
```

The second sum becomes simple.

Since:

```math
s - tk - k = s - (t + 1)k
```

we get:

```math
\sum_{k \ge B} DP[s - tk - k][k]
=
\sum_{k \ge B} DP[s - (t + 1)k][k]
=
x[s][t + 1]
```

This term increases the slope from `t` to `t + 1`.

The first sum shifts the line diagonally.

Let:

```math
j = k - 1
```

Then:

```math
k = j + 1
```

So the term becomes:

```math
DP[s - t(j + 1) - 1][j]
```

which equals:

```math
DP[s - t - 1 - tj][j]
```

For:

```math
j \ge B
```

this contributes:

```math
x[s - t - 1][t]
```

One boundary term appears when:

```math
j = B - 1
```

This belongs to the small-`k` region already computed.

Its contribution is:

```math
DP[s - tB - 1][B - 1]
```

Combining the two parts gives the recurrence:

```math
x[s][t]
=
x[s - t - 1][t]
+
DP[s - tB - 1][B - 1]
+
x[s][t + 1]
```

All invalid indices contribute zero.

For example, if either condition holds:

```math
s - t - 1 < 0
```

```math
s - tB - 1 < 0
```

then the corresponding term is treated as zero.

This recurrence computes a complete grouped line sum rather than a single DP state.

For each fixed slope `t`, process all values of `s` from small to large.

The term:

```math
x[s - t - 1][t]
```

already exists because it uses a smaller value of `s`.

The term:

```math
x[s][t + 1]
```

already exists if you process `t` in decreasing order.

So the evaluation order is:

```text
1. Process t from B down to 0.
2. For each fixed t, process s from small to large.
```

You only need:

```math
t \le B
```

A term in `x[s][t]` contributes only when:

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

Once `t` grows beyond about `sqrt(N)`, the term:

```math
s - tk
```

becomes negative for all useful values of `s` satisfying:

```math
s \le N
```

Such states contribute zero.

So the number of useful slopes is:

```math
O(\sqrt{N})
```

For each slope, you process `O(N)` values of `s`.

Therefore, the grouped large-`k` region costs:

```math
O(N\sqrt{N})
```

The small-`k` region also costs:

```math
O(N\sqrt{N})
```

The total complexity is:

```math
O(N\sqrt{N})
```

The memory stays linear if you keep only the current and next slope layers of `x`, plus the precomputed small-`k` boundary values.

The needed boundary values are:

```math
DP[n][B - 1]
```

The grouped recurrence only needs two layers:

```math
x[\cdot][t]
```

and:

```math
x[\cdot][t + 1]
```

So the memory usage becomes:

```math
O(N)
```

The core idea is simple.

Compute the small columns directly.

For the large columns, do not visit every cell.

Sum cells along structured lines.

Then derive a recurrence for those line sums from the original DP recurrence.

This reduces the computation from:

```math
O(N^2)
```

to:

```math
O(N\sqrt{N})
```
