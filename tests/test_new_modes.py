"""(d) The modes added in November 2026 (ext, tries, mult) against brute force in Python."""
import math

import goldbach_ref as g

D = 20000


def _odd_primes(isp, lim):
    return [a for a in range(3, lim, 2) if isp[a]]


def _covers(n, a, isp):
    return len(g.good_pairs(n, a, isp)) > 0


def _is_p2(x, isp):
    """Omega(x) <= 2 (x >= 2)."""
    c = 0
    d = 2
    while d * d <= x:
        while x % d == 0:
            x //= d
            c += 1
        d += 1
    if x > 1:
        c += 1
    return 1 <= c <= 2


def _chen_covers(n, a, isp):
    m = n - a
    for p in range(3, m, 2):
        if (p & m) != p:
            continue
        q = n - p
        if (isp[p] and _is_p2(q, isp)) or (isp[q] and _is_p2(p, isp)):
            return True
    return False


def _parse_ext(out):
    res = {"crec": [], "upto": {}}
    for line in out.splitlines():
        t = line.split()
        if t and t[0] == "crecord":
            res["crec"].append((int(t[1]), int(t[2])))
        elif t and t[0] == "upto":
            kv = dict(zip(t[2::2], t[3::2]))
            res["upto"][int(t[1])] = kv
        elif t and t[0] == "chen_hard3_upto_1e6":
            res["chen_list"] = [int(x) for x in t[2:]]
    return res


def test_ext_against_python(astar):
    isp = g.sieve(D)
    P = _odd_primes(isp, 1000)
    hard = chenhard = binfail = 0
    chen_list, crec, best = [], [], -1
    for n in range(10, D + 1, 2):
        A = [a for a in P if a + 5 <= n]
        ks = next(i for i, a in enumerate(A) if _covers(n, a, isp))
        kc = next((i for i, a in enumerate(A[:ks]) if _chen_covers(n, a, isp)), ks)
        hard += ks > 0
        if kc > 0:
            chenhard += 1
            chen_list.append(n)
        if kc > best:
            best = kc
            crec.append((n, A[kc]))
        abin = max(A, key=lambda a: (bin(n - a).count("1"), -a))
        binfail += not _covers(n, abin, isp)
    r = _parse_ext(astar("ext", D, 4))
    u = r["upto"][10000]
    # the C program reports cumulative counts at powers of ten; compare at 10^4 and recompute there
    h4 = sum(1 for n in range(10, 10001, 2) if not _covers(n, 3, isp))
    assert int(u["hard3"]) == h4
    assert r["chen_list"] == chen_list
    assert r["crec"] == crec
    assert binfail == 0 and int(r["upto"][D]["binfail"]) == 0
    assert int(r["upto"][D]["hard3"]) == hard and int(r["upto"][D]["chenhard3"]) == chenhard


def test_lemma_low_bits():
    """Lemma: low k bits of n - a are 2^k - e, popcount formula."""
    for n in range(16, 5000, 2):
        for k in range(2, int(math.log2(n)) + 1):
            d = n % (1 << k)
            for a in range(d + 1, d + (1 << k) + 1):
                e = a - d
                assert (n - a) % (1 << k) == (1 << k) - e
                assert bin(n - a).count("1") == bin(n // (1 << k) - 1).count("1") + k - bin(e - 1).count("1")


def test_mult_against_python(astar):
    N = 3000
    isp = g.sieve(N)
    out = astar("mult", "full", 10, N)
    lines = {int(t[1]): t for t in (l.split() for l in out.splitlines()) if t and t[0] == "n"}
    for n in range(10, N + 1, 2):
        pairs = [p for p in range(3, n - 2, 2) if isp[p] and isp[n - p]]
        M = {p: sum(1 for a in range(3, n - p, 2) if isp[a] and ((n - a) & p) == p) for p in pairs}
        t = lines[n]
        assert int(t[2]) == len(pairs)
        assert int(t[3]) == sum(1 for p in pairs if M[p] > 0)
        assert int(t[6]) == sum(M.values())
    # examples quoted in the text
    assert lines[14][2:4] == ["3", "1"]


def test_tries_small(astar):
    out = astar("tries", 200000)
    assert "nofound 0" in out
    assert "mean_tries" in out and " bin 1.0000" in out
