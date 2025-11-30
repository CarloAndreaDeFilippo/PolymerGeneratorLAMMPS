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

        self.sigma_bead    = float(par["sigma_bead"])
        self.sigma_patch   = float(par["sigma_patch"])
        self.sigma_solvent = float(par["sigma_solvent"])
        self.sigma_colloid = float(par["sigma_colloid"])

        if self.nsolvent > 0:
            self.atom_types += 1
            self.solventType = self.atom_types

        if self.npatch > 0:
            self.atom_types += 1
            self.bond_types += 2
            self.patchType = self.atom_types

        if self.ncolloids > 0:
            self.atom_types += 1
            self.colloidType = self.atom_types

        self.max_sigma = max(self.sigma_bead, self.sigma_patch, self.sigma_solvent, self.sigma_colloid)

        self.spheres: list[Sphere] = []*self.ntot
        self.atomList = LinkedCellList(self.Lbox, self.max_sigma, self.ntot)

        self.coordsFirstBeads = []
        self.patchyBeadsIDs = []         #List of tuples AtomID Bead + Patch to keep track of the IDs for the bonds

        self.atomID=1
        self.molID=1

    def addPolymers(self):

        minPolymerDist = self.sigma_bead

        if self.npatch > 0:
            minPolymerDist +=  2. * self.sigma_patch

        distBetweenBeads = 1 * self.sigma_bead #TODO: Define the distance between the beads depending on the chosen potential

        for p in range(self.npol):

            percent = 100 * (p + 1) / self.npol
            print(f"\rPolymer {p+1}/{self.npol}  ({percent:6.2f}%)", end="")

            #First bead of the polymers
            ID_start = int(1 + p * (self.ns + self.npatch))

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
            patchyBeads = sorted(np.random.choice(self.ns, size=self.npatch, replace=False).tolist())

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
                theta = np.random.rand() * 2. * m.pi

                patchCoord[0] += m.cos(theta) * 0.5 * (self.sigma_bead + self.sigma_patch)
                patchCoord[1] += m.sin(theta) * 0.5 * (self.sigma_bead + self.sigma_patch)

                patchCoordPBC = [c - hl for c, hl in zip(patchCoord, self.LboxHalf)]

                for ax in range(3):
                    patchCoordPBC[ax] -= self.Lbox[ax] * round(patchCoordPBC[ax] / self.Lbox[ax])

                self.spheres.append(Sphere(self.sigma_patch, patchCoordPBC, self.atomID, self.molID, self.patchType))

                self.atomList.addObjectToList(self.atomID - 1, [c - hl for c, hl in zip(patchCoordPBC, self.LboxHalf)])

                #Save the tuple AtomID of the bead + patch for the bond

                self.patchyBeadsIDs.append((ID_start + nparticle, self.atomID))

                self.atomID += 1
        
        self.molID += 1
        print()


    #TODO: merge addColloids and addSolvent methods
    def addColloids(self):

        for nc in range(self.ncolloids):

            percent = 100 * (nc + 1) / self.ncolloids
            print(f"\rColloid {nc+1}/{self.ncolloids}  ({percent:6.2f}%)", end="")

            while True:

                colloid = Sphere(self.sigma_colloid, 
                                (np.random.uniform(-self.LboxHalf[0], self.LboxHalf[0]),
                                np.random.uniform(-self.LboxHalf[1], self.LboxHalf[1]),
                                np.random.uniform(-self.LboxHalf[2], self.LboxHalf[2])),
                                self.atomID,
                                self.molID,
                                self.colloidType)

                if self.atomList.overlapCheck(colloid, self.spheres) == False:
                    break

            self.spheres.append(colloid)
            
            self.atomList.addObjectToList(self.atomID - 1, colloid.cm)
            
            self.atomID += 1
            self.molID += 1

        print()

    def addSolvent(self):

        for nc in range(self.nsolvent):

            percent = 100 * (nc + 1) / self.nsolvent
            print(f"\rSolvent {nc+1}/{self.nsolvent}  ({percent:6.2f}%)", end="")

            while True:

                solvent = Sphere(self.sigma_solvent, 
                                (np.random.uniform(-self.LboxHalf[0], self.LboxHalf[0]),
                                np.random.uniform(-self.LboxHalf[1], self.LboxHalf[1]),
                                np.random.uniform(-self.LboxHalf[2], self.LboxHalf[2])),
                                self.atomID,
                                self.molID,
                                self.solventType)

                if self.atomList.overlapCheck(solvent, self.spheres) == False:
                    break
            
            self.spheres.append(solvent)

            self.atomList.addObjectToList(self.atomID - 1, solvent.cm)
            
            self.atomID += 1
            self.molID += 1

        print()
    
    def writeLAMMPSfile(self):

        with open(self.par["atomFile"], "w+") as f:

            f.write("LAMMPS data file\n\n")

            f.write(f"{self.ntot} atoms\n")
            f.write(f"{self.nbonds} bonds\n")
            f.write("0 angles\n0 dihedrals\n0 impropers\n\n")

            f.write(f"{self.atom_types} atom types\n{self.bond_types} bond types\n0 angle types\n0 dihedral types\n0 improper types\n\n")

            f.write(f"0.0 {self.Lbox[0]} xlo xhi\n0.0 {self.Lbox[1]} ylo yhi\n0.0 {self.Lbox[2]} zlo zhi\n\n")

            f.write(f"Masses\n\n{self.beadType} 1\n")

            if self.nsolvent > 0:
                f.write(f"{self.solventType} 0.5\n")  #!Solvent mass 0.5 #TODO: read from par file

            if self.npatch > 0:
                f.write(f"{self.patchType} 1\n")   #!Patch mass 1 #TODO: read from par file

            if self.ncolloids > 0:
                f.write(f"{self.colloidType} 1\n")     #!Colloid mass 1 #TODO: read from par file

            f.write("\n")

            f.write("Atoms\n")
            f.write("# atom-ID\tmol-ID\tatom-type\tx\ty\tz\n")

            for sph in self.spheres:
                f.write(f"{sph.atomID} {sph.molID} {sph.atomType} {sph.cm[0]+self.LboxHalf[0]} {sph.cm[1]+self.LboxHalf[1]} {sph.cm[2]+self.LboxHalf[2]}\n")

            f.write("\n")


            # -----------> VELOCITIES <------------- # #*It is possible to specify the atoms velocities

            #f.write("Velocities\n\n")

            #[f.write("{0} {1} {2} {3}\n".format(i, np.random.random() * 2., np.random.random()  * 2., np.random.random()  * 2. - 1.)) for i in range(1, ntot + 1)]

            #f.write("\n")


            # -----------> BONDS <------------- #

            # Bond IDs
            f.write("Bonds\n#bond-ID\tbond-type\tfirst-atom\tsecond-atom\n")

            ID_bond = 1

            #Bonds between the beads
            for p in range(self.npol):

                #First bead of the polymer
                ID_start = int(1 + p * (self.ns + self.npatch))

                for i in range(self.ns - 1):

                    atom_id = ID_start + i

                    f.write(f"{ID_bond} 1 {atom_id} {atom_id + 1}\n")

                    ID_bond += 1

            #Bonds between the bead and the patch
            for BeadID,PatchID in self.patchyBeadsIDs:

                f.write(f"{ID_bond} 2 {BeadID} {PatchID}\n")
                ID_bond += 1


    def saveXYZfile(self) -> None:

        namesVMD = {1: 'C', 2:'H', 3:'O', 4:'N'}

        with open(self.par["xyzFile"], "w+") as f:
            f.write(f"{self.ntot}\n\n")
            for sphere in self.spheres:
                letter = namesVMD[sphere.atomType]
                f.write(f"{letter} {sphere.cm[0]} {sphere.cm[1]} {sphere.cm[2]}\n")