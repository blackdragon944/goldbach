"""(c) Reproduce the previously known numbers."""
import goldbach_ref as g
from conftest import parse_scan

KNOWN_RECORDS = [(10, 3), (12, 5), (38, 7), (302, 19), (556, 47), (2618, 59), (2642, 83),
                 (16502, 127), (786638, 211), (8388826, 223), (16843016, 269), (201326986, 397)]


def test_records_and_no_exceptions_up_to_2_2e8(astar):
    r = parse_scan(astar("scan", 220_000_000, 4))
    assert r["exceptions"] == [8]
    assert r["records"] == KNOWN_RECORDS


def test_histogram_1e5(astar):
    r = parse_scan(astar("scan", 100_000, 4))
    h = r["hist"][100_000]
    assert (h[3], h[5], h[7], h[11]) == (48058, 1037, 349, 280)


def test_histogram_1e5_python():
    t = g.a_star_table(100_000)
    vals = [a for n, a in t.items() if n >= 10]
    assert [vals.count(a) for a in (3, 5, 7, 11)] == [48058, 1037, 349, 280]
    assert t[8] is None and t[6] is None


def test_family_2k_plus_4_fails_for_a3():
    isp = g.sieve(1 << 21)
    for k in range(3, 21):
        assert not g.has_good_pair((1 << k) + 4, 3, isp)


def test_carries_variant_2e6(astar):
    out = astar("carries", 2_000_000, 4)
    got = {int(l.split()[0][len("carries<="):]): int(l.split()[2]) for l in out.splitlines()}
    assert got == {0: 14284, 1: 2446, 2: 476, 3: 110}


def test_68_exception_for_all_c():
    isp = g.sieve(100)
    assert g.min_carries_a3(68, isp) > 3


def test_even_a_useless():
    """n - a even, p odd => lowest bit of p not in n - a."""
    for n in range(10, 200, 2):
        for a in range(2, n, 2):
            for p in range(3, n, 2):
                assert (p & (n - a)) != p
