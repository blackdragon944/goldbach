"""
Reference (slow, simple) implementation of the "Sierpinski-Goldbach" condition.

Condition for even n, odd prime a:
    there exist odd primes p, q with p + q = n, a < q, and C(n - a, p) odd.

Equivalent forms (see tests/test_equivalences.py):
    Lucas   : (p & (n - a)) == p          ("p is a binary submask of n - a")
    Kummer  : p + (q - a) has no carries in base 2
    Geometry: cell (p, 2^K + a - 1 - n) of the Sierpinski triangle T_a is filled

a*(n) = least odd prime a for which the condition holds (None if no such a).

This file is deliberately written for clarity, not speed.  The fast program
is src/astar.c; tests/test_c_vs_python.py checks that the two agree.
"""

from __future__ import annotations

import math
from functools import lru_cache


# ---------------------------------------------------------------- primes

def sieve(limit: int) -> bytearray:
    """is_prime[i] for 0 <= i <= limit (bytearray of 0/1)."""
    s = bytearray([1]) * (limit + 1)
    s[0:2] = b"\x00\x00"[: min(2, limit + 1)]
    for i in range(2, math.isqrt(limit) + 1):
        if s[i]:
            s[i * i :: i] = bytearray(len(range(i * i, limit + 1, i)))
    return s


def odd_primes_upto(limit: int, isp: bytearray) -> list[int]:
    return [i for i in range(3, limit + 1, 2) if isp[i]]


# ---------------------------------------------------------- the condition

def is_submask(p: int, m: int) -> bool:
    return (p & m) == p


def carries(x: int, y: int) -> int:
    """Number of carries when adding x + y in base 2 (Kummer: = v_2(C(x+y, x)))."""
    c = 0
    count = 0
    while x or y or c:
        s = (x & 1) + (y & 1) + c
        c = 1 if s >= 2 else 0
        count += c
        x >>= 1
        y >>= 1
    return count


def good_pairs(n: int, a: int, isp: bytearray) -> list[int]:
    """All odd primes p with q = n - p an odd prime, q > a, and p a submask of n - a.

    (q > a is the same as p < n - a; together with p ⊆ n - a it means p != n - a.)
    """
    m = n - a
    out = []
    for p in range(3, m, 2):
        if (p & m) == p and isp[p] and isp[n - p]:
            out.append(p)
    return out


def has_good_pair(n: int, a: int, isp: bytearray) -> bool:
    """Same as bool(good_pairs(...)) but enumerates submasks of m = n - a directly."""
    m = n - a
    if m <= 0:
        return False
    p = (m - 1) & m  # largest proper submask of m
    while p > 0:
        if p & 1 and p >= 3 and isp[p] and isp[n - p]:
            return True
        p = (p - 1) & m
    return False


def a_star(n: int, isp: bytearray, primes: list[int]) -> int | None:
    """Least odd prime a with a good pair, or None (only n = 8 among even n >= 6 so far)."""
    for a in primes:
        if a >= n - 3:  # need q > a and q <= n - 3
            break
        if has_good_pair(n, a, isp):
            return a
    return None


def a_star_table(N: int) -> dict[int, int | None]:
    """a*(n) for all even n with 6 <= n <= N."""
    isp = sieve(N)
    primes = odd_primes_upto(N, isp)
    return {n: a_star(n, isp, primes) for n in range(6, N + 1, 2)}


# ------------------------------------------------- variant: allow c carries

def min_carries_a3(n: int, isp: bytearray) -> int | None:
    """min over odd primes p, q = n - p prime, q > 3, of carries(p, q - 3)."""
    best = None
    for p in range(3, n - 4, 2):
        q = n - p
        if q > 3 and isp[p] and isp[q]:
            c = carries(p, q - 3)
            if best is None or c < best:
                best = c
                if c == 0:
                    break
    return best


# ---------------------------------------------- geometry of the triangle T_a

def filled_cells(a: int, K: int):
    """Filled cells (u, v) of T_a:  (u & v) == 0 and u + v <= 2^K - 1."""
    S = 1 << K
    for u in range(S):
        for v in range(S - u):
            if (u & v) == 0:
                yield u, v


def hypotenuses(a: int, K: int) -> set[tuple[tuple[int, int], tuple[int, int]]]:
    """Hypotenuse segments of all filled half-cells of T_a, in plane coordinates.

    Cell (u, v) occupies x in [u, u+1], y in [y0, y0+1] with y0 = 2^K + a - 1 - v;
    its filled upper-left half is the triangle (u, y0), (u, y0+1), (u+1, y0+1),
    whose hypotenuse runs from (u, y0) to (u+1, y0+1) along a line of slope 1.
    """
    S = 1 << K
    segs = set()
    for u, v in filled_cells(a, K):
        y0 = S + a - 1 - v
        segs.add(((u, y0), (u + 1, y0 + 1)))
    return segs


def geometric_condition(n: int, a: int, p: int, segs) -> bool:
    """Does the line y = x + q (q = n - p) run along a side of a filled triangle at (p, n)?"""
    return ((p, n), (p + 1, n + 1)) in segs


if __name__ == "__main__":
    import sys

    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    t = a_star_table(N)
    for n, a in t.items():
        print(n, a)
