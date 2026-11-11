import os
import subprocess
import sys

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
EXE = os.path.join(ROOT, "bin", "astar.exe" if os.name == "nt" else "astar")


def run_astar(*args) -> str:
    if not os.path.exists(EXE):
        pytest.skip("bin/astar not built (see README)")
    out = subprocess.run([EXE, *map(str, args)], capture_output=True, text=True, check=True)
    return out.stdout


def parse_scan(text: str) -> dict:
    res = {"records": [], "hist": {}, "exceptions": []}
    for line in text.splitlines():
        if line.startswith("record "):
            _, n, a = line.split()
            res["records"].append((int(n), int(a)))
        elif line.startswith("exceptions"):
            res["exceptions"] = [int(x) for x in line.split(":", 1)[1].split()]
        elif line.startswith("hist n<="):
            head, rest = line.split(":", 1)
            upto = int(head[len("hist n<="):])
            h = {}
            for tok in rest.split():
                k, v = tok.split(":")
                h[k if k == "big" else int(k)] = int(v)
            res["hist"][upto] = h
    return res


@pytest.fixture(scope="session")
def astar():
    return run_astar
