import numpy as np
import math as m
import re
from os import listdir, path

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

                self.atomList.addObjectToList(self.atomID - 1, patchCoordPBC)

                #Save the tuple AtomID of the bead + patch for the bond

                self.patchyBeadsIDs.append((ID_start + nparticle, self.atomID))

                self.atomID += 1
        
        self.molID += 1
        print()


    def addSpheres(self, nSpheres: int, sigma: float, sphereType: int, sphereName: str):

        for nc in range(nSpheres):

            percent = 100 * (nc + 1) / nSpheres
            print(f"\r{sphereName} {nc+1}/{nSpheres}  ({percent:6.2f}%)", end="")

            while True:

                sphere = Sphere(sigma, 
                                (np.random.uniform(-self.LboxHalf[0], self.LboxHalf[0]),
                                np.random.uniform(-self.LboxHalf[1], self.LboxHalf[1]),
                                np.random.uniform(-self.LboxHalf[2], self.LboxHalf[2])),
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

    
    def writeConfigurationFileLAMMPS(self):

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


    @staticmethod
    def minPositionLJ(sigma: float) -> float:
        """Minimum of LJ to mimick WCA"""
        return 2**(1.0/6.0) * sigma

    @staticmethod
    def maxDistanceInteractionLJ(sigma: float) -> float:
        """Maxmimum distance interaction for LJ"""
        return 2.5 * sigma

    def writeInputFileLAMMPS(self):

        with open(f"{self.par["inputFile"]}", "w+") as f:

            f.write("######################\n#   Initialization   #\n######################\n\n")
	  
            f.write("units lj\n")
            f.write("atom_style bond\n")
            f.write("special_bonds fene\n") #To not count the LJ potential twice for bonded beads
            f.write("boundary p p p\n")
            f.write("comm_modify mode single cutoff 5.0 vel yes")   #For thread communication (mpi)

            f.write("\n")

            f.write("#######################\n#   Atom definition   #\n#######################\n\n")


            if int(self.par['restart']) == 1:
                restart_files = [f for f in listdir("restarts") if f.endswith(".restart")]
                restart_paths = [path.join("restarts", f) for f in restart_files]
                latestRestart = max(restart_paths, key=path.getsize)

                f.write(f"read_restart {latestRestart} remap\n")
                f.write("\n")
            else:
                f.write(f"read_data {self.par["atomFile"]} extra/special/per/atom 100\n") # da modificare per prendere parametri esterni
                f.write("\n")

            f.write("################\n#   Settings   #\n################\n\n")

            #*Spheres setup

            #Cut and shifted LJ potential at the minimum (mimic WCA potential)
            f.write(f"pair_style lj/cut {self.minPositionLJ(self.par['sigma_bead'])}\n")
            f.write(f"pair_modify shift yes\n")

            #LJ interaction: type1 type2 epsilon sigma cutoff

            #Beads-beads interaction
            f.write(f"pair_coeff {self.beadType} {self.beadType} {self.par['eps']} {self.par['sigma_bead']} {self.minPositionLJ(self.par['sigma_bead'])}\n")   #! WCA

            if self.par['nsolvent'] > 0:
                #Bead-solvent interaction
                sigma_bead_solv = 0.5 * (self.par['sigma_bead'] + self.par['sigma_solvent'])
                f.write(f"pair_coeff {self.beadType} {self.solventType} {self.par['eps']} {sigma_bead_solv} {self.minPositionLJ(sigma_bead_solv)}\n")   #! WCA
                
                #Solvent-solvent interaction
                f.write(f"pair_coeff {self.solventType} {self.solventType} {self.par['eps_ss']} {self.par['sigma_solvent']} {self.maxDistanceInteractionLJ(self.par['sigma_solvent'])}\n") #! LJ cut

            if self.par['npatch'] > 0:
                #Bead-patch interaction
                sigma_bead_patch = 0.5 * (self.par['sigma_bead'] + self.par['sigma_patch'])
                f.write(f"pair_coeff {self.beadType} {self.patchType} {self.par['eps']} {sigma_bead_patch} {self.minPositionLJ(sigma_bead_patch)}\n")   #! WCA
                
                #Patch-patch interaction
                f.write(f"pair_coeff {self.patchType} {self.patchType} {self.par['eps']} {self.par['sigma_patch']} {self.minPositionLJ(self.par['sigma_patch'])}\n")   #! WCA

                if self.par['nsolvent'] > 0:
                    #Patch-solvent
                    sigma_patch_solv = 0.5 * (self.par['sigma_patch'] + self.par['sigma_solvent'])
                    f.write(f"pair_coeff {self.patchType} {self.solventType} {self.par['eps']} {sigma_patch_solv} {self.minPositionLJ(sigma_patch_solv)}\n")   #! WCA


            if self.par['ncolloids'] > 0:
                #Bead-colloid interaction
                sigma_bead_coll = 0.5 * (self.par['sigma_bead'] + self.par['sigma_colloid'])
                f.write(f"pair_coeff {self.beadType} {self.colloidType} {self.par['eps_bc']} {sigma_bead_coll} {self.maxDistanceInteractionLJ(sigma_bead_coll)}\n")   #! LJ cut
                
                #Colloid-colloid interaction
                f.write(f"pair_coeff {self.colloidType} {self.colloidType} {self.par['eps']} {self.par['sigma_colloid']} {self.minPositionLJ(self.par['sigma_colloid'])}\n")    #! WCA

                if self.par['nsolvent'] > 0:
                    #Colloid-solvent interaction
                    sigma_coll_solv = 0.5 * (self.par['sigma_colloid'] + self.par['sigma_solvent'])
                    f.write(f"pair_coeff {self.colloidType} {self.solventType} {self.par['eps_cs']} {sigma_coll_solv} {self.maxDistanceInteractionLJ(sigma_coll_solv)}\n")   #! LJ cut
            
                if self.par['npatch'] > 0:
                    #Colloid-patch interaction
                    sigma_coll_patch = 0.5 * (self.par['sigma_colloid'] + self.par['sigma_patch'])
                    f.write(f"pair_coeff {self.colloidType} {self.patchType} {self.par['eps']} {sigma_coll_patch} {self.minPositionLJ(sigma_coll_patch)}\n")   #! WCA
                

            f.write("\n")

            #*Bonds setup

            f.write("neigh_modify one 10000\n")

            f.write("bond_style fene\n")

            #FENE bond: bondID Kparameter MaxDistance LJeps sigma

            #Bond between beads
            f.write(f"bond_coeff 1 30.0 {1.5 * self.par['sigma_bead']} {self.par['eps']} {self.par['sigma_bead']}\n") 

            if self.par['npatch'] > 0:
                #Bond between bead and patch
                sigma_bead_patch = 0.5 * (self.par['sigma_bead'] + self.par['sigma_patch'])
                f.write(f"bond_coeff 2 30.0 {1.5 * sigma_bead_patch} {self.par['eps']} {sigma_bead_patch}\n")

                #Bond for the crosslink between patches
                f.write(f"bond_coeff 3 30.0 {1.5 * self.par['sigma_patch']} {self.par['eps']} {self.par['sigma_patch']}\n")

            f.write("\n")	 


            #*Simulation environment setup

            f.write("fix 1 all nve\n")                             #NVE integrator (combined with Langevin thermostat gives NVT)
            f.write("fix 2 all langevin 1.0 1.0 1.0 16113\n")      #Langevin thermostat (Tstart Tfin dampening seed)

            if self.par['press/berendsen'] > 0:
                f.write(f"fix 3 all press/berendsen iso {self.par['press_in']} {self.par['press_fin']} 1000.0\n") #If NPT with Berendsen barostat

            #f.write("velocity all create 50.0 66254 rot yes dist gaussian\n") #Random initial velocities



            f.write("#####################\n#   Calculations    #\n#####################\n")
            f.write("\n")

            file_number = 0

            if int(self.par['restart']) == 1:
                    
                    lammpsTrjFiles = [f for f in listdir("configurations") if "polymers" in f]

                    numbers = [int(re.search(r"(\d+)", f).group(1)) for f in lammpsTrjFiles]

                    file_number = max(numbers, default=-1) + 1


            #TODO: fix the computation of gyration radius for the polymers

            '''
            f.write("compute molchunk all chunk/atom molecule\n") #Definisco i chunk in base al molID

            f.write("compute rg all gyration/chunk molchunk\n")
            f.write("compute com all com/chunk molchunk\n")

            f.write("fix save_rg all ave/time 100 1 100 c_rg file gyration/gyration.txt mode vector\n")
            f.write("fix save_com all ave/time 100 1 100 c_com[*] file com/com.txt mode vector\n")
            '''

            #f.write("compute gyrb all gyration\n")
            #f.write(f"fix gyration all ave/time 1000 1 1000 c_gyrb file gyration/gyr_{file_number}.txt start 0\n\n")

            #f.write("compute comb all com\n")
            #f.write(f"fix com all ave/time 1000 1 1000 c_comb[*] file com/com_{file_number}.txt start 0\n\n")

            '''
            f.write("##Group polymers\n")

            for p in range(0, 1):

                #primo atomo del polimero
                ID_start = int(1 + p * (par['ns'] + par['npatch']))
                ID_end = int(1 + (p + 1) * (par['ns'] + par['npatch']))

                f.write(f"group polymer{p} id {ID_start}:{ID_end}\n")
                f.write(f"compute gyrb{p} polymer{p} gyration\n")
                f.write(f"fix gyration polymer{p} ave/time 1000 1 1000 c_gyrb{p} file gyration/gyr_{file_number}_{p}.txt start 0\n\n")
                
                #f.write(f"compute com{p} polymer{p} com\n")
                #f.write(f"fix com polymer{p} ave/time 1000 50 500000 c_com{p} file com/com_{file_number}_{p}.txt start 500000\n\n")
            '''
            #!compute com0 polymer0 com

            '''
            # Rg backbone
            f.write("group bbone id {0}:{1}\n".format(np.min(bbone_ids), np.max(bbone_ids)))
            f.write("compute gyrb bbone gyration\n")
            f.write("fix gyration bbone ave/time 1000 50 500000 c_gyrb file gyration/gyr_bbone_{}.txt start 500000\n".format(file_number)) # da modficare per prendere parametri esterni
            f.write("\n")
            '''

            # Rg complete
            #f.write("compute gyr all gyration\n")
            #f.write("fix gyration1 all ave/time 1000 50 50000 c_gyr file gyration/gyr_all_{}.txt start 50000\n".format(file_number))
            #f.write("\n")

            '''
            # End-to-end
            f.write("group bbone_extremes id {0} {1}\n".format(np.min(bbone_ids), np.max(bbone_ids)))
            f.write("\n")
            '''


            f.write("###########\n# Dumping #\n###########\n")
            f.write("\n")

            f.write(f"restart {int(self.par["restartSteps"])} restarts/polymers.restart restarts/polymers2.restart\n")  #LAMMPS restart files
            f.write("\n")
            f.write(f"dump config all custom {int(self.par['dumpsteps'])} configurations/polymers_{file_number}.lammpstrj id mol type x y z\n")
            f.write("\n")
            #f.write("dump endtoend bbone_extremes xyz 500000 gyration/bbone_extremes_*.lammpstrj\n")# id type x y z\n".format(file_number))
            #f.write("\n")

            f.write("##########################\n#   Running simulation   #\n##########################\n")
            f.write("\n")
            f.write("timestep 0.001\n") #Timestep of the integrator 
            f.write(f"thermo {int(self.par['thermosteps'])}\n")
            f.write("thermo_style custom step temp etotal ke pe epair emol press vol\n")
            #!thermo_style custom step temp etotal ke pe epair emol press vol c_com0[1] c_com0[2] c_com0[3]
            f.write("thermo_modify flush yes\n")

            #*Crosslink: bond creation during the simulation between patches
            #fix fixName groupID bond/create Nevery itype jtype Rmin bondtype prob p seed iparam maxbonds newitype jparam maxbonds newjtype

            #Create the crosslinks if there are patches in the system and the crosslink parameter is set to True (value != 0)
            if self.par['npatch'] > 0 and self.par['crosslink'] == True:
                f.write(f"fix createBonds all bond/create 1 {self.patchType} {self.patchType} {1 * self.par['sigma_bead']} 3 iparam 1 2 jparam 1 2\n")

            f.write("minimize 1.0e-4 1.0e-6 100 1000\n") #Little energy minimization to stabilize initial configuration

            f.write(f"run {int(self.par['totsteps'])}\n")
