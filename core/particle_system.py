import numpy as np
import math as m

from .sphere import Sphere
from .linked_cell_list import LinkedCellList

class ParticleSystem:

    def __init__(self, par: dict):

        self.par = par

        self.Lbox = par["Lbox"]
        self.LboxHalf = [0.5*L for L in self.Lbox]

        self.npol = int(par["npol"])
        self.ns = int(par["ns"])
        self.npatch = int(par["npatch"])
        self.nsolvent = int(par["nsolvent"])
        self.ncolloids = int(par["ncolloids"])

        self.ntot = self.npol*(self.ns + self.npatch) + self.nsolvent + self.ncolloids

        self.nbonds = self.npol * (self.ns - 1 + self.npatch)

        self.atom_types = 1
        self.bond_types = 1

        self.beadType = 1

        self.sigma_bead    = float(par["bead"]["sigma"])
        self.sigma_patch   = float(par["patch"]["sigma"])
        self.sigma_solvent = float(par["solvent"]["sigma"])
        self.sigma_colloid = float(par["colloid"]["sigma"])

        sigmas = [self.sigma_bead]

        if self.nsolvent > 0:
            self.atom_types += 1
            self.solventType = self.atom_types
            sigmas.append(self.sigma_solvent)

        if self.npatch > 0:
            self.atom_types += 1
            self.bond_types += 2
            self.patchType = self.atom_types
            sigmas.append(self.sigma_patch)

        if self.ncolloids > 0:
            self.atom_types += 1
            self.colloidType = self.atom_types
            sigmas.append(self.sigma_colloid)

        self.max_sigma = max(sigmas)

        self.spheres: list[Sphere] = []
        self.atomList = LinkedCellList(self.Lbox, self.max_sigma, self.ntot)

        self.coordsFirstBeads = []
        self.patchyBeadsIDs = []         #List of tuples AtomID Bead + Patch to keep track of the IDs for the bonds

        self.atomID=1
        self.molID=1

        self.rng = np.random.default_rng()

    def addPolymers(self):

        minPolymerDist = self.sigma_bead

        if self.npatch > 0:
            minPolymerDist +=  2. * self.sigma_patch

        distBetweenBeads = 1 * self.sigma_bead #TODO: Define the distance between the beads depending on the chosen potential

        for p in range(self.npol):

            print(f"\rPolymer {p+1}/{self.npol}  ({100 * (p + 1) / self.npol:6.2f}%)", end="")

            #First bead of the polymers
            ID_start = int(1 + p * (self.ns + self.npatch))

            #Random coordinate of first bead (not overlapping with other polymers)
            while True:

                coordFirstBead = [self.rng.random() * l for l in self.Lbox]

                for coord in self.coordsFirstBeads:

                    dist = [abs(ci - cj) for (ci, cj) in zip(coord, coordFirstBead)]

                    for ax in range(2):

                        if dist[ax] >= self.LboxHalf[ax]:
                            dist[ax] = dist[ax] - self.Lbox[ax]
                    
                    distance2DSqrd = dist[0]*dist[0] + dist[1]*dist[1]

                    if distance2DSqrd <= minPolymerDist*minPolymerDist:
                        break
                        
                else:   #all the beads are placed without overlap
                    break

            self.coordsFirstBeads.append(coordFirstBead)

            #List of the indexes of patch beads (sorted to have the correct molID)
            patchyBeads = sorted(self.rng.choice(self.ns, size=self.npatch, replace=False).tolist())

            #Place the beads
            for npart in range(self.ns):

                coord = coordFirstBead.copy()
                coord[2] += distBetweenBeads * npart   

                coordPBC = [c - hl for c, hl in zip(coord, self.LboxHalf)]

                for ax in range(3):
                    coordPBC[ax] -= self.Lbox[ax] * round(coordPBC[ax] / self.Lbox[ax])

                self.spheres.append(Sphere(self.sigma_bead, coordPBC, self.atomID, self.molID, self.beadType))
                self.atomList.addObjectToList(self.atomID - 1, coordPBC)

                self.atomID += 1
    
            #Place the patches

            for nparticle in patchyBeads:
                
                patchCoord = coordFirstBead.copy()
                patchCoord[2] += distBetweenBeads * nparticle

                #Random position of the patch around the bead
                theta = self.rng.random() * 2. * m.pi

                patchCoord[0] += m.cos(theta) * 0.5 * (self.sigma_bead + self.sigma_patch)
                patchCoord[1] += m.sin(theta) * 0.5 * (self.sigma_bead + self.sigma_patch)

                patchCoordPBC = [c - hl for c, hl in zip(patchCoord, self.LboxHalf)]

                for ax in range(3):
                    patchCoordPBC[ax] -= self.Lbox[ax] * round(patchCoordPBC[ax] / self.Lbox[ax])

                self.spheres.append(Sphere(self.sigma_patch, patchCoordPBC, self.atomID, self.molID, self.patchType))

                self.atomList.addObjectToList(self.atomID - 1, patchCoordPBC)

                #Save the tuple AtomID of the bead + patch for the bond

                self.patchyBeadsIDs.append((ID_start + nparticle, self.atomID))

                self.atomID += 1
        
        self.molID += 1
        print()


    def addSpheres(self, nSpheres: int, sigma: float, sphereType: int, sphereName: str):

        for nc in range(nSpheres):

            print(f"\r{sphereName} {nc+1}/{nSpheres}  ({100 * (nc + 1) / nSpheres:6.2f}%)", end="")

            while True:

                sphere = Sphere(sigma, 
                                (self.rng.uniform(-self.LboxHalf[0], self.LboxHalf[0]),
                                self.rng.uniform(-self.LboxHalf[1], self.LboxHalf[1]),
                                self.rng.uniform(-self.LboxHalf[2], self.LboxHalf[2])),
                                self.atomID,
                                self.molID,
                                sphereType)

                if self.atomList.overlapCheck(sphere, self.spheres) == False:
                    break

            self.spheres.append(sphere)
            
            self.atomList.addObjectToList(self.atomID - 1, sphere.cm)
            
            self.atomID += 1
            self.molID += 1

        print()

    def addColloids(self):

        self.addSpheres(self.ncolloids, self.sigma_colloid, self.colloidType, "Colloid")

    def addSolvent(self):

        self.addSpheres(self.nsolvent, self.sigma_solvent, self.solventType, "Solvent")
    