## The \(O(N\sqrt{N})\) DP Trick for Partition Numbers

This trick applies to a DP with states \(DP[n][k]\), where:

\[
DP[n][k]
\]

means the number of ways to form \(n\) using exactly \(k\) positive parts.

The standard recurrence is:

\[
DP[n][k] = DP[n-1][k-1] + DP[n-k][k]
\]

The term \(DP[n-1][k-1]\) handles the case where each part loses \(1\).

The term \(DP[n-k][k]\) keeps the number of parts equal to \(k\), and removes \(1\) from each of the \(k\) parts.

A direct implementation fills all states \(DP[n][k]\). This costs:

\[
O(N^2)
\]

We reduce this to:

\[
O(N\sqrt{N})
\]

---

## Split the DP Table

Choose:

\[
B = \lfloor \sqrt{N} \rfloor
\]

Split the DP table into two regions.

The first region is:

\[
k < B
\]

The second region is:

\[
k \ge B
\]

For \(k < B\), compute the DP directly.

This region has \(B\) columns and \(N\) rows, so it costs:

\[
O(NB) = O(N\sqrt{N})
\]

Keep these values. The large-\(k\) region will need the boundary column:

\[
DP[n][B-1]
\]

---

## Group the Large-\(k\) States

For \(k \ge B\), avoid computing each state alone.

Instead, group states along lines.

Define:

\[
x[s][t] = \sum_{k \ge B} DP[s - tk][k]
\]

Here:

- \(s\) is the starting value.
- \(t\) is the slope parameter.
- \(k\) is the DP column.
- \(DP[s - tk][k]\) is one point on the line.

The large-\(k\) contribution to the answer is:

\[
x[N][0]
\]

because:

\[
x[N][0] = \sum_{k \ge B} DP[N][k]
\]

So the final partition number is:

\[
p(N) =
\sum_{1 \le k < B} DP[N][k]
+
x[N][0]
\]

---

## Why Only \(O(\sqrt{N})\) Slopes Matter

A term inside \(x[s][t]\) is valid only when:

\[
s - tk \ge 0
\]

Since \(k \ge B\), we have:

\[
tk \ge tB
\]

When \(t\) becomes larger than about \(\sqrt{N}\), the value \(s - tk\) becomes negative for all useful states.

So we only need:

\[
0 \le t \le B
\]

This gives \(O(\sqrt{N})\) slopes.

---

## Derive the Recurrence for \(x[s][t]\)

Start from:

\[
x[s][t] = \sum_{k \ge B} DP[s - tk][k]
\]

Use the original recurrence:

\[
DP[n][k] = DP[n-1][k-1] + DP[n-k][k]
\]

Apply it to each term:

\[
x[s][t]
=
\sum_{k \ge B} DP[s - tk - 1][k - 1]
+
\sum_{k \ge B} DP[s - tk - k][k]
\]

Now simplify both sums.

The second sum becomes:

\[
\sum_{k \ge B} DP[s - (t+1)k][k]
=
x[s][t+1]
\]

For the first sum, substitute:

\[
j = k - 1
\]

Then:

\[
DP[s - tk - 1][k - 1]
=
DP[s - t(j+1) - 1][j]
\]

So:

\[
DP[s - t(j+1) - 1][j]
=
DP[s - t - 1 - tj][j]
\]

For \(j \ge B\), these terms form:

\[
x[s - t - 1][t]
\]

There is one boundary term when:

\[
j = B - 1
\]

This term belongs to the small-\(k\) region already computed. It contributes:

\[
DP[s - tB - 1][B - 1]
\]

So the grouped recurrence is:

\[
x[s][t]
=
x[s - t - 1][t]
+
DP[s - tB - 1][B - 1]
+
x[s][t+1]
\]

Invalid indices contribute zero.

---

## Evaluation Order

Compute \(x[s][t]\) with this order:

1. Process \(t\) from \(B\) down to \(0\).
2. For each fixed \(t\), process \(s\) from small to large.

This works because:

\[
x[s][t]
\]

depends on:

\[
x[s - t - 1][t]
\]

which has a smaller \(s\), and:

\[
x[s][t+1]
\]

which belongs to the next slope layer.

So you only need to keep two slope layers in memory.

---

## Complexity

The small-\(k\) region costs:

\[
O(N\sqrt{N})
\]

The grouped large-\(k\) region has:

\[
O(\sqrt{N})
\]

slope values.

For each slope, we process:

\[
O(N)
\]

values of \(s\).

So the large-\(k\) region also costs:

\[
O(N\sqrt{N})
\]

The total time complexity is:

\[
O(N\sqrt{N})
\]

The memory usage is:

\[
O(N)
\]

if we keep only the needed DP columns and two layers of \(x\).

---

## Core Idea

Compute small \(k\) directly.

For large \(k\), group DP states along lines:

\[
DP[s - tk][k]
\]

Then derive a recurrence for these line sums using the original DP recurrence.

This avoids visiting all large-\(k\) cells one by one.
