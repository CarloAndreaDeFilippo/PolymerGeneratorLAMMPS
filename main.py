import argparse

from generate_polymers import generate_polymers
from generate_input import generate_input

from core.parse_parameters import parseParameters
from core.particle_system import ParticleSystem

def main():
    
    parser = argparse.ArgumentParser(prog="generate_input.py")
    parser.add_argument(
        "parameter_file",
        help="Path to the parameter file (json)."
    )

    args = parser.parse_args()
    parFile = args.parameter_file

    partsys = ParticleSystem(parseParameters(parFile))

    generate_polymers(partsys)
    generate_input(partsys)

if __name__ == "__main__":
    main()
