#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

from core.parse_parameters import parseParameters
from core.particle_system import ParticleSystem

def generate_polymers(partsys: ParticleSystem):

    print("Generating configuration file...")
    
    partsys.addPolymers()
    partsys.addColloids()
    partsys.addSolvent()

    partsys.writeConfigurationFileLAMMPS()

    if partsys.par["saveXYZ"] == True:
        partsys.saveXYZfile()

def main():

    parser = argparse.ArgumentParser(prog="generate_polymers.py")
    parser.add_argument(
        "parameter_file",
        help="Path to the parameter file (json)."
    )

    args = parser.parse_args()
    parFile = args.parameter_file

    partsys = ParticleSystem(parseParameters(parFile))

    generate_polymers(partsys)


if __name__ == "__main__":
    main()