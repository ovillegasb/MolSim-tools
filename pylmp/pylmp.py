"""
Conjunto de funciones para trabajar con la topologia, la construccion de sistemas
y la exportacion de la informacion.

"""

from pylmp.inputBOSS import genzmat


def save_topol(file, res):
    """Save topology."""

    # save RES.z
    genzmat(file, res)
