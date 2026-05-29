# The $O(N\sqrt{N})$ DP Trick for Partition Numbers

Use this optimization when your DP table has two coordinates, $n$ and $k$, and the transition moves in two regular directions.

For the partition DP, define $DP[n][k]$ as the number of partitions of $n$ into exactly $k$ positive parts.

The standard recurrence is:

$$
DP[n][k] = DP[n-1][k-1] + DP[n-k][k]
$$

The first term, $DP[n-1][k-1]$, reduces each part by $1$.

The second term, $DP[n-k][k]$, keeps the number of parts equal to $k$, and removes $1$ from each of the $k$ parts.

A direct solution fills all states $DP[n][k]$ for $1 \le n \le N$ and $1 \le k \le N$. This costs:

$$
O(N^2)
$$

The optimization starts by splitting the DP table at:

$$
B = \lfloor \sqrt{N} \rfloor
$$

For small values of $k$, meaning $k < B$, compute the DP normally.

This region has $B$ columns and $N$ rows, so its cost is:

$$
O(NB) = O(N\sqrt{N})
$$

This small-$k$ region gives all values $DP[n][k]$ for $k < B$.

Keep these values because the large-$k$ computation will touch their boundary.

For large values of $k$, meaning $k \ge B$, avoid computing every state separately.

Group many states together along lines.

Define:

$$
x[s][t] = \sum_{k \ge B} DP[s - tk][k]
$$

Here, $t$ represents the slope of the line.

For fixed $s$ and $t$, the expression sums all large-$k$ DP values lying on one line in the DP table.

The value needed from the large-$k$ region is:

$$
x[N][0] = \sum_{k \ge B} DP[N][k]
$$

So the final answer becomes:

$$
p(N) =
\sum_{1 \le k < B} DP[N][k]
+
x[N][0]
$$

Now apply the original recurrence to every term inside $x[s][t]$.

Start with:

$$
x[s][t] =
\sum_{k \ge B} DP[s - tk][k]
$$

Using:

$$
DP[n][k] = DP[n-1][k-1] + DP[n-k][k]
$$

we get:

$$
x[s][t]
=
\sum_{k \ge B} DP[s - tk - 1][k - 1]
+
\sum_{k \ge B} DP[s - tk - k][k]
$$

The second sum becomes simple:

$$
\sum_{k \ge B} DP[s - (t+1)k][k]
=
x[s][t+1]
$$

The first sum shifts the line diagonally.

Let:

$$
j = k - 1
$$

Then the term becomes:

$$
DP[s - t(j+1) - 1][j]
$$

which equals:

$$
DP[s - t - 1 - tj][j]
$$

For $j \ge B$, this contributes:

$$
x[s - t - 1][t]
$$

One boundary term appears when $j = B - 1$.

This belongs to the small-$k$ region already computed.

Its contribution is:

$$
DP[s - tB - 1][B - 1]
$$

So the recurrence for the grouped values is:

$$
x[s][t]
=
x[s - t - 1][t]
+
DP[s - tB - 1][B - 1]
+
x[s][t+1]
$$

Treat invalid indices as zero.

This recurrence turns many large-$k$ DP states into one grouped transition.

For each slope $t$, compute all $s$ from small to large.

The term $x[s - t - 1][t]$ already exists because it uses a smaller $s$.

The term $x[s][t+1]$ already exists if you process $t$ in decreasing order.

You only need $t \le B$.

Since $k \ge B$, the term $s - tk$ becomes negative once $t$ grows beyond about $\sqrt{N}$.

Such states contribute zero.

So the number of useful slopes is:

$$
O(\sqrt{N})
$$

For each slope, you process $O(N)$ values of $s$.

Therefore, the grouped large-$k$ region costs:

$$
O(N\sqrt{N})
$$

The total complexity is:

$$
O(N\sqrt{N})
$$

The memory stays linear if you keep only the current and next slope layers of $x$, plus the precomputed small-$k$ boundary values.

The core idea is simple.

Compute the small columns directly.

For the large columns, do not visit every cell.

Sum cells along structured lines, then derive a recurrence for those line sums from the original DP recurrence.
