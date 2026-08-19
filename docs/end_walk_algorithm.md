# `end_walk` algorithm (v0.2.2 audit)

**Coordinates:** BED 0-based half-open `[start, end)`.  
**Series:** one numeric value per bin.  
**Jump array:** `d[k] = |y[k+1] - y[k]|` lives on the *edge between*
bin `k` and bin `k+1` (`k = 1, …, n-1`).

## Why v0.2.1 was directional

The v0.2.1 implementation (`version = "legacy"`):

```text
search d[i0 - r, i0 + r] around EACH end independently
return argmax index k
L = start[k]
R = end[k]
if R-index < L-index: swap
```

Three asymmetries:

1. **Same mapping for both ends.** A jump at edge `k` is the change
   between bins `k` and `k+1`. The high block that *starts* after that
   jump has left coordinate `start[k+1]`, not `start[k]`. The block that
   *ends* at that jump has right coordinate `end[k]`. Using `start[k]`
   for `L` shifts the left end one bin upstream.
2. **`which.max` on a wide window.** Search radius was the whole
   extract half-window (often `2 ×` domain width), so both the left and
   right jumps were visible from either end. `which.max` returns the
   *first* maximum, i.e. the left jump, for both searches.
3. **Swap instead of per-end assignment.** If the right search returns
   a smaller index than the left search, the two ends are swapped. That
   is not a symmetric left/right procedure.

These are implementation issues, not evidence that Gaussian change-in-mean
is intrinsically left-biased.

## v0.2.2 symmetric procedure

Search radius `r` is local to each candidate end and is **not** allowed
to cross the midpoint of the candidate.

### Left search

```text
mid = floor((iL + iR) / 2)
edges = {k : max(1, iL - r) ≤ k ≤ min(mid, iL + r)}
kL = argmax d[k] on edges
      (ties: rightmost, i.e. toward the candidate interior)
L  = start[kL + 1]     # first bin to the right of the jump
```

### Right search

```text
edges = {k : max(mid, iR - r) ≤ k ≤ min(n-1, iR + r)}
kR = argmax d[k] on edges
      (ties: leftmost, i.e. toward the candidate interior)
R  = end[kR]           # last bin to the left of the jump
```

If `R ≤ L`, fall back to the candidate `[start[iL], end[iR])`.

The two searches are mirror images: reverse the series, reverse the
interval as `[N-R, N-L)`, and the mapped errors swap.

## What this is not

- Not a peak caller.
- Not an NB likelihood.
- Not a claim that broad, low-gradient domains will have zero error
  under noise. Broad failure after the oracle test is a local-method
  limitation if the oracle candidate still has large error.
