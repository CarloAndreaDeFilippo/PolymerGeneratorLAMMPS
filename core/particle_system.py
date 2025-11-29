import numpy as np
import math as m

from .sphere import Sphere
from .linked_cell_list import LinkedCellList

class ParticleSystem:

    coordsFirstBeads = []
    patchyBeadsIDs = []         #List of tuples AtomID Bead + Patch to keep track of the IDs for the bonds

    atomID=1
    molID=1

    def __init__(self, Lbox: list[float], max_sigma: float, ntot: int):
        self.Lbox = Lbox
        self.LboxHalf = [0.5*L for L in Lbox]
        self.spheres: list[Sphere] = []*ntot
        self.atomList = LinkedCellList(Lbox, max_sigma, ntot)

    def addPolymers(self, npol:int, ns:int, npatch:int, sigma_bead:float, sigma_patch:float, beadType:int, patchType:int):

        minPolymerDist = sigma_bead

        if npatch > 0:
            minPolymerDist +=  2. * sigma_patch

        distBetweenBeads = 1 * sigma_bead #TODO: Define the distance between the beads depending on the chosen potential

        for p in range(npol):

            percent = 100 * (p + 1) / npol
            print(f"\rPolymer {p+1}/{npol}  ({percent:6.2f}%)", end="")

            #First bead of the polymers
            ID_start = int(1 + p * (ns + npatch))

            #Random coordinate of first bead (not overlapping with other polymers)
            while True:

                coordFirstBead = [np.random.random() * l for l in self.Lbox]

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
            patchyBeads = sorted(np.random.choice(ns, size=npatch, replace=False).tolist())

            #Place the beads
            for npart in range(ns):

                coord = coordFirstBead.copy()
                coord[2] += distBetweenBeads * npart   

                coordPBC = [c - hl for c, hl in zip(coord, self.LboxHalf)]

                for ax in range(3):
                    coordPBC[ax] -= self.Lbox[ax] * round(coordPBC[ax] / self.Lbox[ax])

                self.spheres.append(Sphere(sigma_bead, coordPBC, self.atomID, self.molID, beadType))
                self.atomList.addObjectToList(self.atomID - 1, coordPBC)

                self.atomID += 1
    
            #Place the patches

            for nparticle in patchyBeads:
                
                patchCoord = coordFirstBead.copy()
                patchCoord[2] += distBetweenBeads * nparticle

                #Random position of the patch around the bead
                theta = np.random.rand() * 2. * m.pi

                patchCoord[0] += m.cos(theta) * 0.5 * (sigma_bead + sigma_patch)
                patchCoord[1] += m.sin(theta) * 0.5 * (sigma_bead + sigma_patch)

                patchCoordPBC = [c - hl for c, hl in zip(patchCoord, self.LboxHalf)]

                for ax in range(3):
                    patchCoordPBC[ax] -= self.Lbox[ax] * round(patchCoordPBC[ax] / self.Lbox[ax])

                self.spheres.append(Sphere(sigma_patch, patchCoordPBC, self.atomID, self.molID, patchType))

                self.atomList.addObjectToList(self.atomID - 1, [c - hl for c, hl in zip(patchCoordPBC, self.LboxHalf)])

                #Save the tuple AtomID of the bead + patch for the bond

                self.patchyBeadsIDs.append((ID_start + nparticle, self.atomID))

                self.atomID += 1
        
            print()
        
        self.molID += 1


    #TODO: merge the two methods
    def addColloids(self, ncoll:int, sigma_colloid:float, colloidType:int):

        for nc in range(ncoll):

            percent = 100 * (nc + 1) / ncoll
            print(f"\rColloid {nc+1}/{ncoll}  ({percent:6.2f}%)", end="")

            while True:

                colloid = Sphere(sigma_colloid, 
                                (np.random.uniform(-self.LboxHalf[0], self.LboxHalf[0]),
                                np.random.uniform(-self.LboxHalf[1], self.LboxHalf[1]),
                                np.random.uniform(-self.LboxHalf[2], self.LboxHalf[2])),
                                self.atomID,
                                self.molID,
                                colloidType)

                if self.atomList.overlapCheck(colloid, self.spheres) == False:
                    break

            self.spheres.append(colloid)
            
            self.atomList.addObjectToList(self.atomID - 1, colloid.cm)
            
            self.atomID += 1
            self.molID += 1

        print()

    def addSolvent(self, nsolv:int, sigma_solvent:float, solventType:int):

        for nc in range(nsolv):

            percent = 100 * (nc + 1) / nsolv
            print(f"\rSolvent {nc+1}/{nsolv}  ({percent:6.2f}%)", end="")

            while True:

                solvent = Sphere(sigma_solvent, 
                                (np.random.uniform(-self.LboxHalf[0], self.LboxHalf[0]),
                                np.random.uniform(-self.LboxHalf[1], self.LboxHalf[1]),
                                np.random.uniform(-self.LboxHalf[2], self.LboxHalf[2])),
                                self.atomID,
                                self.molID,
                                solventType)

                if self.atomList.overlapCheck(solvent, self.spheres) == False:
                    break
            
            self.spheres.append(solvent)

            self.atomList.addObjectToList(self.atomID - 1, solvent.cm)
            
            self.atomID += 1
            self.molID += 1

        print()