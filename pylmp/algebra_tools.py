"""
Module for compiling mathematical functions.

"""

import numpy as np
from numpy import linalg as lg


def Distance(u, v):
    """Return length of a vector."""

    u = np.array(u)
    v = np.array(v)

    return lg.norm(u - v)


def pairing_func(a, b):
    ans = (a + b) * (a + b + 1) * 0.5
    if a > b:
        ans = ans + a
    else:
        ans = ans + b
    return (int(ans))
