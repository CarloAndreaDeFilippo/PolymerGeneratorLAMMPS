from .particle_system import ParticleSystem

def writeLAMMPSConfiguration(partsys: ParticleSystem):
    with open(partsys.par["atomFile"], "w+") as f:

            f.write("LAMMPS data file\n\n")

            f.write(f"{partsys.ntot} atoms\n")
            f.write(f"{partsys.nbonds} bonds\n")
            f.write("0 angles\n0 dihedrals\n0 impropers\n\n")

            f.write(f"{partsys.atom_types} atom types\n{partsys.bond_types} bond types\n0 angle types\n0 dihedral types\n0 improper types\n\n")

            f.write(f"0.0 {partsys.Lbox[0]} xlo xhi\n0.0 {partsys.Lbox[1]} ylo yhi\n0.0 {partsys.Lbox[2]} zlo zhi\n\n")

            f.write(f"Masses\n\n{partsys.beadType} {float(partsys.par["bead"]["mass"])}\n")

            if partsys.nsolvent > 0:
                f.write(f"{partsys.solventType} {float(partsys.par["solvent"]["mass"])}\n")

            if partsys.npatch > 0:
                f.write(f"{partsys.patchType} {float(partsys.par["patch"]["mass"])}\n")

            if partsys.ncolloids > 0:
                f.write(f"{partsys.colloidType} {float(partsys.par["colloid"]["mass"])}\n")

            f.write("\n")

            f.write("Atoms\n")
            f.write("# atom-ID\tmol-ID\tatom-type\tx\ty\tz\n")

            for sph in partsys.spheres:
                f.write(f"{sph.atomID} {sph.molID} {sph.atomType} {sph.cm[0]+partsys.LboxHalf[0]} {sph.cm[1]+partsys.LboxHalf[1]} {sph.cm[2]+partsys.LboxHalf[2]}\n")

            f.write("\n")


            # -----------> VELOCITIES <------------- # #*It is possible to specify the atoms velocities

            if partsys.par["randomVelocities"]["enabled"] == True:

                v = partsys.par["randomVelocities"]["maxMagnitude"]

                f.write("Velocities\n\n")

                [f.write(f"{i} {partsys.rng.random() * v[0]} {partsys.rng.random() * v[1]} {partsys.rng.random() * v[2]}\n") for i in range(1, partsys.ntot + 1)]

                f.write("\n")


            # -----------> BONDS <------------- #

            # Bond IDs
            f.write("Bonds\n#bond-ID\tbond-type\tfirst-atom\tsecond-atom\n")

            ID_bond = 1

            #Bonds between the beads
            for p in range(partsys.npol):

                #First bead of the polymer
                ID_start = int(1 + p * (partsys.ns + partsys.npatch))

                for i in range(partsys.ns - 1):

                    atom_id = ID_start + i

                    f.write(f"{ID_bond} 1 {atom_id} {atom_id + 1}\n")

                    ID_bond += 1

            #Bonds between the bead and the patch
            for BeadID,PatchID in partsys.patchyBeadsIDs:

                f.write(f"{ID_bond} 2 {BeadID} {PatchID}\n")
                ID_bond += 1
