# An O(N sqrt(N)) Dynamic Programming Optimization for Partition Numbers

This section presents an optimization for a class of two-dimensional dynamic programming recurrences indexed by two parameters, \(n\) and \(k\). The method applies when the transition structure follows two regular directions in the DP table. The partition-number recurrence provides a clear example of this setting.

Let \(DP[n][k]\) denote the number of partitions of \(n\) into exactly \(k\) positive parts. The standard recurrence is given by:

```math
DP[n][k] = DP[n - 1][k - 1] + DP[n - k][k]
```

The first term, \(DP[n - 1][k - 1]\), corresponds to reducing each part by \(1\). This removes one unit from every part and decreases both the sum and the number of parts by one. The second term, \(DP[n - k][k]\), keeps the number of parts fixed at \(k\), while removing one unit from each of the \(k\) parts. Together, these two cases form the standard recurrence for partitions into exactly \(k\) positive parts.

A direct computation evaluates every state \(DP[n][k]\) for \(1 \le n \le N\) and \(1 \le k \le N\). This requires:

```math
O(N^2)
```

time. The optimization reduces this complexity by splitting the DP table into two regions according to the value of \(k\). Let:

```math
B = \lfloor \sqrt{N} \rfloor
```

The first region consists of states with \(k < B\). This region has \(B\) columns and \(N\) rows. It is small enough to compute directly using the original recurrence. The cost of this step is:

```math
O(NB) = O(N\sqrt{N})
```

The values computed in this small-\(k\) region must be retained, since the large-\(k\) computation reaches the boundary column \(DP[n][B - 1]\). This boundary appears when a transition from a state with column \(k\) moves to column \(k - 1\), and crosses from the large-\(k\) region into the precomputed small-\(k\) region.

The second region consists of states with \(k \ge B\). Computing this region state by state would still be expensive. Instead, the method groups many large-\(k\) states along structured lines in the DP table. Define the auxiliary quantity:

```math
x[s][t] = \sum_{k \ge B} DP[s - tk][k]
```

Here, \(s\) represents the starting coordinate, while \(t\) represents the slope parameter. For fixed \(s\) and \(t\), the expression sums all large-\(k\) DP states lying on one line of the form \(DP[s - tk][k]\). This aggregation replaces many individual states with one grouped value.

The contribution of the large-\(k\) region to the final partition number is obtained from:

```math
x[N][0] = \sum_{k \ge B} DP[N][k]
```

Therefore, the partition number \(p(N)\) is expressed as the sum of the directly computed small-\(k\) contribution and the grouped large-\(k\) contribution:

```math
p(N)
=
\sum_{1 \le k < B} DP[N][k]
+
x[N][0]
```

It remains to derive an efficient recurrence for \(x[s][t]\). Starting from the definition:

```math
x[s][t]
=
\sum_{k \ge B} DP[s - tk][k]
```

and applying the original recurrence to each term gives:

```math
x[s][t]
=
\sum_{k \ge B} DP[s - tk - 1][k - 1]
+
\sum_{k \ge B} DP[s - tk - k][k]
```

The second summation has a direct interpretation. Since \(s - tk - k = s - (t + 1)k\), it follows that:

```math
\sum_{k \ge B} DP[s - tk - k][k]
=
\sum_{k \ge B} DP[s - (t + 1)k][k]
=
x[s][t + 1]
```

Thus, this part of the recurrence moves from slope \(t\) to slope \(t + 1\).

The first summation shifts the line diagonally. Let \(j = k - 1\). Then \(k = j + 1\), and the summand becomes:

```math
DP[s - t(j + 1) - 1][j]
```

Equivalently:

```math
DP[s - t - 1 - tj][j]
```

For all \(j \ge B\), these terms form the grouped value:

```math
x[s - t - 1][t]
```

There is one boundary term when \(j = B - 1\). This term lies outside the large-\(k\) region and belongs to the precomputed small-\(k\) region. Its contribution is:

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

All invalid indices contribute zero. For example, if \(s - t - 1 < 0\) or \(s - tB - 1 < 0\), the corresponding term is treated as zero.

The recurrence above computes a complete grouped line sum rather than a single DP state. For each fixed slope \(t\), the values of \(s\) are processed in increasing order. This ensures that the dependency \(x[s - t - 1][t]\) has already been computed. The slopes are processed in decreasing order of \(t\), which ensures that \(x[s][t + 1]\) is already available.

The number of useful slope values is bounded by \(O(\sqrt{N})\). A term in \(x[s][t]\) contributes only if \(s - tk \ge 0\). Since the grouped region satisfies \(k \ge B\), we have \(tk \ge tB\). Once \(t\) exceeds approximately \(\sqrt{N}\), the expression \(s - tk\) becomes negative for all relevant \(s \le N\). Therefore, only slopes \(0 \le t \le B\) need to be processed.

The small-\(k\) region requires \(O(NB)\) time. The grouped large-\(k\) region processes \(O(B)\) slopes, and each slope requires \(O(N)\) work. Its cost is also \(O(NB)\). Since \(B = \lfloor \sqrt{N} \rfloor\), the total time complexity is:

```math
O(N\sqrt{N})
```

The memory usage remains linear when the implementation stores only the required DP columns, the boundary values \(DP[n][B - 1]\), and two consecutive slope layers of \(x\). Thus, the method achieves:

```math
O(N)
```

memory.

The main idea is to compute the small-\(k\) columns directly and replace the large-\(k\) cells with line sums. The original recurrence then induces a recurrence over these line sums. This removes the need to visit every large-\(k\) state individually and reduces the computation from \(O(N^2)\) to \(O(N\sqrt{N})\).
