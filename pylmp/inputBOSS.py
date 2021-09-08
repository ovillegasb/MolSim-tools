"""
Modulo dedicado a prepara el archivo de entrada para ser tratado por BOSS.

"""

import os
import pylmp.algebra_tools as alg
import networkx as nx


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


def make_graphs(atoms, coos, bonds):
    G = nx.DiGraph()

    # Add nodes using atom types and coordinate.
    for i in coos.keys():
        G.add_node(
            i,
            XYZ=coos[i],
            elem=atoms[i],
            atno=alg.bossElement2Num(atoms[i])
        )

    for (i, j, rij) in zip(bonds['BI'], bonds['BJ'], bonds['RIJ']):
        G.add_edge(i, j, distance=rij)
        G.add_edge(j, i, distance=rij)

    all_ps = dict(nx.algorithms.all_pairs_shortest_path_length(G))
    all_paths = []

    for s in all_ps.keys():
        for e in all_ps[s].keys():
            if all_ps[s][e] == 1:
                all_paths += list(nx.algorithms.all_simple_paths(G, s, e, cutoff=1))

            elif all_ps[s][e] == 2:
                all_paths += list(nx.algorithms.all_simple_paths(G, s, e, cutoff=2))

            elif all_ps[s][e] == 3:
                all_paths += list(nx.algorithms.all_simple_paths(G, s, e, cutoff=3))

    all_bonds = [p for p in all_paths if len(set(p)) == 2]
    new_angs = [p for p in all_paths if len(set(p)) == 3]
    new_tors = [p for p in all_paths if len(set(p)) == 4]

    dict_new_tors = {alg.tor_id(t): t for t in new_tors}
    dict_new_angs = {alg.ang_id(t): t for t in new_angs}

    imp_keys = [n for n in G.nodes() if G.degree(n) / 2 == 3]
    all_imps = {}
    for i in imp_keys:
        nei = list(G.neighbors(i))
        if G.nodes[i]['atno'] == 6:
            all_imps[i] = [nei[0], i, nei[1], nei[2]]
    MOL_ICOORDS = {'BONDS': all_bonds,
                   'ANGLES': dict_new_angs,
                   'TORSIONS': dict_new_tors,
                   'IMPROPERS': all_imps
                   }

    print(MOL_ICOORDS)
    exit()

    return G, MOL_ICOORDS


def get_zmat(file, res):
    """Build a initial file RES.z from input file.mol."""

    # read mol lines
    lines = open(file, 'r').readlines()

    # extract information from mol file
    coord, atypes, bonds = readmolfile(lines)

    # Extract connectivity using networkx
    G_mol, mol_icords = make_graphs(atypes, coord, bonds)

    print('RES.z saved!')


def genzmat(file, res):

    # Log file (olog)
    os.system("rm -vf olog")
    os.system("touch olog")

    # save RES.z initial
    get_zmat(file, res)
    pass
