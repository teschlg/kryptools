"""
Discrete log solvers: Shanks' baby-step-giant-step algorithm
"""

from math import isqrt


def dlog_bsgs(a: int, b: int, n: int, m: int = None) -> int:
    """Compute the discrete log_a(b) in Z_n of an element a of order m using Shanks' baby-step-giant-step algorithm."""
    a %= n
    b %= n
    if not m:
        m = n - 1
    mm = 1 + isqrt(m - 1)
    # initialize baby_steps table
    baby_steps = {}
    baby_step = 1
    for k in range(mm):
        baby_steps[baby_step] = k
        baby_step = baby_step * a % n

    # now take the giant steps
    giant_stride = pow(a, -mm, n)
    giant_step = b
    for l in range(mm):
        if giant_step in baby_steps:
            return l * mm + baby_steps[giant_step]
        giant_step = giant_step * giant_stride % n
    raise ValueError("DLP not solvable.")
