"""(b) The C program and the Python reference give the same a*(n) for all even n <= 10^6."""
import os
import tempfile

import numpy as np

import goldbach_ref as g
from conftest import parse_scan

N = 10**6


def test_astar_c_equals_python(astar):
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "astar.bin")
        astar("scan", N, 4, N, f)
        c = np.fromfile(f, dtype=np.uint16)
    py = g.a_star_table(N)
    bad = [n for n, a in py.items() if n >= 8 and int(c[n // 2]) != (a or 0)]
    assert not bad, bad[:20]
    assert py[8] is None and c[4] == 0


def test_G3_c_equals_python(astar):
    """G_3(n) (number of good pairs) agrees for n <= 3000."""
    D = 3000
    isp = g.sieve(D)
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "g3.bin")
        astar("count", 3, D, f)
        G = np.fromfile(f, dtype=np.uint32)
    for n in range(8, D + 1, 2):
        assert G[n // 2] == len(g.good_pairs(n, 3, isp)), n


def test_carries_c_equals_python(astar):
    """Carry-variant exception counts agree with brute force for n <= 3000."""
    D = 3000
    isp = g.sieve(D)
    mins = [g.min_carries_a3(n, isp) for n in range(10, D + 1, 2)]
    out = astar("carries", D, 4)
    for line in out.splitlines():
        c = int(line.split()[0][len("carries<="):])
        e = int(line.split()[2])
        assert e == sum(1 for x in mins if x is None or x > c), (c, e)
