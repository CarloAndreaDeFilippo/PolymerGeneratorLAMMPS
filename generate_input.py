#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from os import getcwd
import argparse

from core.parse_parameters import parseParameters
from core.make_directories import makeDirectories

from core.particle_system import ParticleSystem

def generate_input(partsys: ParticleSystem):

    print("Generating input file...")

    #Creation of the necessary directories
    makeDirectories([getcwd() + "/configurations", getcwd() + "/restarts"])

    partsys.writeInputFileLAMMPS()


def main():

    parser = argparse.ArgumentParser(prog="generate_input.py")
    parser.add_argument(
        "parameter_file",
        help="Path to the parameter file (json)."
    )

    args = parser.parse_args()
    parFile = args.parameter_file

    partsys = ParticleSystem(parseParameters(parFile))

    generate_input(partsys)

if __name__ == "__main__":
    main()