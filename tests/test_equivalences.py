"""(a) Geometric, bitwise (Lucas), binomial-parity and Kummer forms agree."""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import goldbach_ref as g  # noqa: E402


def test_four_definitions_agree():
    isp = g.sieve(200)
    odd_primes = g.odd_primes_upto(200, isp)
    checked = 0
    for a in [3, 5, 7, 11, 13]:
        for K in (7, 8):  # 2^K > n - a for all n tested; result must not depend on K
            segs = g.hypotenuses(a, K)
            for n in range(a + 5, 120, 2):  # a odd, so a + 5 is even
                for p in odd_primes:
                    if p >= n - a:
                        break
                    bitwise = (p & (n - a)) == p
                    parity = math.comb(n - a, p) % 2 == 1
                    kummer = g.carries(p, (n - p) - a) == 0
                    geom = g.geometric_condition(n, a, p, segs)
                    assert bitwise == parity == kummer == geom, (a, K, n, p)
                    checked += 1
    assert checked > 5000


def test_all_integers_not_only_primes():
    """The equivalence is purely combinatorial: test every 0 <= p <= m < 64."""
    K = 6
    for a in (3, 5):
        segs = g.hypotenuses(a, K)
        for m in range(0, 64):
            n = m + a
            for p in range(0, m + 1):
                bitwise = (p & m) == p
                parity = math.comb(m, p) % 2 == 1
                kummer = g.carries(p, m - p) == 0
                geom = g.geometric_condition(n, a, p, segs)
                assert bitwise == parity == kummer == geom, (a, m, p)


def test_trivial_when_a_equals_q():
    """a = q: p = n - a is always a submask of itself."""
    for n in range(10, 400, 2):
        for q in (3, 5, 7, 11):
            assert ((n - q) & (n - q)) == n - q


def test_n_14_figure_pairs():
    """Pairs for the figure: 3+11 on a side (a=3), 7+7 in a hole, 11+3 trivial (a = q)."""
    n, a = 14, 3
    m = n - a  # 11 = 1011b
    assert (3 & m) == 3          # 3 + 11: good
    assert (7 & m) != 7          # 7 + 7 : hole
    assert (11 & m) == 11        # 11 + 3: q = a, excluded
    segs = g.hypotenuses(a, 4)
    assert g.geometric_condition(n, a, 3, segs)
    assert not g.geometric_condition(n, a, 7, segs)
