from .sphere import Sphere

def saveXYZ(fileName : str, spheres : list[Sphere]) -> None:

    namesVMD = {1: 'C', 2:'H', 3:'O', 4:'N'}

    with open(fileName, "w+") as f:
        f.write(f"{len(spheres)}\n\n")
        for sphere in spheres:
            letter = namesVMD[sphere.atomType]
            f.write(f"{letter} {sphere.cm[0]} {sphere.cm[1]} {sphere.cm[2]}\n")