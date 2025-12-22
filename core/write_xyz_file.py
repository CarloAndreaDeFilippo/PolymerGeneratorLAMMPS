
from .particle_system import ParticleSystem

def writeXYZfile(partsys: ParticleSystem) -> None:

    namesVMD = {1: 'C', 2:'H', 3:'O', 4:'N'}

    with open(partsys.par["xyzFile"], "w+") as f:
        f.write(f"{partsys.ntot}\n\n")
        for sphere in partsys.spheres:
            letter = namesVMD[sphere.atomType]
            f.write(f"{letter} {sphere.cm[0]} {sphere.cm[1]} {sphere.cm[2]}\n")
