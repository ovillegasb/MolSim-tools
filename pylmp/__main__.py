#!/usr/bin/env python
# -*- conding=utf-8 -*-

r"""
  _______     ___      __  __ _____
 |  __ \ \   / / |    |  \/  |  __ \
 | |__) \ \_/ /| |    | \  / | |__) |
 |  ___/ \   / | |    | |\/| |  ___/
 | |      | |  | |____| |  | | |
 |_|      |_|  |______|_|  |_|_|

PYLMP parte de la idea del servidor LigParGen para generar estructuras y topologias
moleculares empleando el campo de fuerzas OPLS-AA pero desde lineas de comandos.

Esto nos ayuda a construir un objeto con la informacion de la topologia y estructura
y usando packmol se construye la configuracion inicial del sistema para simulacion
de dinamica molecular en LAMMPS.

Es posible usar openbabel, solo lo usan para convertir formatos

export BOSSdir="/home/ragnar/applications/LigParGen/boss" Asi se exporta a las
variables de entorno.

Existen scripts en boss que nos ayudan a calcular las cargas y a optimizar la molecula.

xZCM1A

"""

import argparse
import time


def options():
    """Generate command line interface"""

    parser = argparse.ArgumentParser(
        prog="PYLMP",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage="%(prog)s [options] files",
        epilog="Enjoy the program!",
        description=__doc__
    )

    # Add the arguments
    # File
    parser.add_argument("-f",
                        help="molecule file coordinate.",
                        default="",
                        metavar="mol",
                        type=str)

    # Name residue
    parser.add_argument("-r",
                        help="Residue name in three words.",
                        default="",
                        metavar="RES",
                        type=str)

    return vars(parser.parse_args())


t0 = time.time()

# starting setups
args = options()


resname = args['r']
file = args['f']

t1 = time.time() - t0
print("Done in %.0f s" % t1)