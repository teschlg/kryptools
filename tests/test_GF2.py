import pytest
from random import randint, seed
from kryptools import GF2
seed(0)

num_tests = 50

def test_GF2_ops():
    for n, modulus in [(4, None), (8, None), (8, 0b100011011), (12, None), (128, 0b11100001 << 120)]:
        gf = GF2(n, modulus=modulus)
        for _ in range(num_tests):
            a = randint(0, gf.order-1)
            b = randint(0, gf.order-1)
            assert (gf(a) + gf(b)).poly() == gf(a).poly() + gf(b).poly()
            assert (gf(a) - gf(b)).poly() == gf(a).poly() - gf(b).poly()
            assert (gf(a) * gf(b)).poly() == gf(a).poly() * gf(b).poly()
            if not (b):
                continue
            assert (gf(a) / gf(b)).poly() == gf(a).poly() / gf(b).poly()
        with pytest.raises(ValueError):
            gf(1) / gf(0)
        with pytest.raises(ValueError):
            gf(0)**-1
        for _ in range(num_tests):
            a = randint(1, gf.order-1)
            assert (gf(a)**-1).poly() == gf(a).poly().inv()
