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


def bossElement2Num(elem):
    symb2mass = {
        'H':  1,
        'B':  5,
        'C':  6,
        'N':  7,
        'O':  8,
        'F':  9,
        'Si': 14,
        'P': 15,
        'S': 16,
        'Cl': 17,
        'Br': 35,
        'I': 53,
    }
    try:
        res = symb2mass[elem]
    except NameError:
        print("Mass for atom %s is not available \n add it to symb2mass dictionary")
    return res


def tor_id(a):
    bond = pairing_func(a[1], a[2])
    ends = pairing_func(a[0], a[3])

    return '%d-%d' % (bond, ends)


def ang_id(a):
    bond_a = pairing_func(a[0], a[1])
    bond_b = pairing_func(a[1], a[2])

    return pairing_func(bond_a, bond_b)
