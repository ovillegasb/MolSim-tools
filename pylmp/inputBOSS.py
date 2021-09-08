"""
Modulo dedicado a prepara el archivo de entrada para ser tratado por BOSS.

"""

import os
import pylmp.algebra_tools as alg


def readmolfile(lines):

    # Search the n atoms and bonds.
    [nats, nbonds] = map(int, (lines[3][0:3], lines[3][3:6]))

    # information lines for atoms
    atlines = lines[4:4 + nats]

    # extract coordinates like dict for n atoms
    # coord[iat] = [xi, yi, zi]
    coord = {}
    atypes = {}
    for i in range(nats):
        els = atlines[i].split()
        coord[i + 1] = [float(e) for e in els[0:3]]
        # atom symbol
        atypes[i + 1] = els[3]

    bondlines = lines[4 + nats:4 + nats + nbonds]
    bonds = {'BI': [], 'BJ': [], 'RIJ': [], 'UID': []}
    for line in bondlines:
        [bi, bj] = map(int, [line[0:3], line[3:6]])
        bonds['BI'].append(bi)
        bonds['BJ'].append(bj)
        bonds['RIJ'].append(alg.Distance(coord[bi], coord[bj]))
        bonds['UID'].append(alg.pairing_func(bi, bj))

    return coord, atypes, bonds


def get_zmat(file, res):
    """Build a initial file RES.z from input file.mol."""

    # read mol lines
    lines = open(file, 'r').readlines()

    # extract information from mol file
    coord, atypes, bonds = readmolfile(lines)

    print('RES.z saved!')


def genzmat(file, res):

    # Log file (olog)
    os.system("rm -vf olog")
    os.system("touch olog")

    # save RES.z initial
    get_zmat(file, res)
    pass
