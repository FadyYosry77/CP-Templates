# The O(N sqrt(N)) DP Trick for Partition Numbers

Use this optimization when your DP table has two coordinates, `n` and `k`, and the transition moves in two regular directions.

For the partition DP, define `DP[n][k]` as the number of partitions of `n` into exactly `k` positive parts.

The standard recurrence is:

```text
DP[n][k] = DP[n - 1][k - 1] + DP[n - k][k]
```

The first term, `DP[n - 1][k - 1]`, reduces each part by `1`.

The second term, `DP[n - k][k]`, keeps the number of parts equal to `k`, and removes `1` from each of the `k` parts.

A direct solution fills all states `DP[n][k]` for `1 <= n <= N` and `1 <= k <= N`.

This costs:

```text
O(N^2)
```

The optimization starts by splitting the DP table at:

```text
B = floor(sqrt(N))
```

For small values of `k`, meaning `k < B`, compute the DP normally.

This region has `B` columns and `N` rows, so its cost is:

```text
O(NB) = O(N sqrt(N))
```

This small-`k` region gives all values `DP[n][k]` for `k < B`.

Keep these values because the large-`k` computation touches their boundary.

For large values of `k`, meaning `k >= B`, avoid computing every state separately.

Group many states together along lines.

Define:

```text
x[s][t] = sum over k >= B of DP[s - t*k][k]
```

Here, `t` represents the slope of the line.

For fixed `s` and `t`, this expression sums all large-`k` DP values lying on one line in the DP table.

The value needed from the large-`k` region is:

```text
x[N][0] = sum over k >= B of DP[N][k]
```

So the final answer becomes:

```text
p(N) = sum over 1 <= k < B of DP[N][k] + x[N][0]
```

Now apply the original recurrence to every term inside `x[s][t]`.

Start with:

```text
x[s][t] = sum over k >= B of DP[s - t*k][k]
```

Using:

```text
DP[n][k] = DP[n - 1][k - 1] + DP[n - k][k]
```

we get:

```text
x[s][t]
=
sum over k >= B of DP[s - t*k - 1][k - 1]
+
sum over k >= B of DP[s - t*k - k][k]
```

The second sum becomes simple:

```text
sum over k >= B of DP[s - (t + 1)*k][k]
=
x[s][t + 1]
```

The first sum shifts the line diagonally.

Let:

```text
j = k - 1
```

Then the term becomes:

```text
DP[s - t*(j + 1) - 1][j]
```

which equals:

```text
DP[s - t - 1 - t*j][j]
```

For `j >= B`, this contributes:

```text
x[s - t - 1][t]
```

One boundary term appears when `j = B - 1`.

This belongs to the small-`k` region already computed.

Its contribution is:

```text
DP[s - t*B - 1][B - 1]
```

So the recurrence for the grouped values is:

```text
x[s][t]
=
x[s - t - 1][t]
+
DP[s - t*B - 1][B - 1]
+
x[s][t + 1]
```

Treat invalid indices as zero.

This recurrence turns many large-`k` DP states into one grouped transition.

For each slope `t`, compute all `s` from small to large.

The term `x[s - t - 1][t]` already exists because it uses a smaller `s`.

The term `x[s][t + 1]` already exists if you process `t` in decreasing order.

You only need `t <= B`.

Since `k >= B`, the term `s - t*k` becomes negative once `t` grows beyond about `sqrt(N)`.

Such states contribute zero.

So the number of useful slopes is:

```text
O(sqrt(N))
```

For each slope, you process `O(N)` values of `s`.

Therefore, the grouped large-`k` region costs:

```text
O(N sqrt(N))
```

The total complexity is:

```text
O(N sqrt(N))
```

The memory stays linear if you keep only the current and next slope layers of `x`, plus the precomputed small-`k` boundary values.

The core idea is simple.

Compute the small columns directly.

For the large columns, do not visit every cell.

Sum cells along structured lines, then derive a recurrence for those line sums from the original DP recurrence.
