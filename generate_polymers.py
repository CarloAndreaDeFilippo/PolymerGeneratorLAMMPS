#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import math as m

from core.sphere import Sphere, SphereType
from core.linked_cell_list import LinkedCellList

# Reading par.txt
par = {}
par_file = 'par.txt'

with open(par_file) as f:
    for line in f:
        (key, val) = line.split()
        par[key] = float(val)

Lbox = [par['lbox_x'], par['lbox_y'], par['lbox_z']]

ntot = int(par['npol'] * (par['ns'] + par['npatch']) + par['nsolvent'])
nbonds = int(par['npol'] * (par['ns'] - 1 + par['npatch']))

if par['ncolloids'] > 0:

    ntot += int(par['ncolloids'])

# ids of all monomers
tot_ids = list(range(1, ntot + 1))

spheres = []*ntot    #Container of all spherical particles in the system

# split of tot_ids into single polymers
#polymer_ids = [tot_ids[x : x + ns] for x in range(0, len(tot_ids), ns)]


f = open("polymers.dat", "w+")

f.write("LAMMPS data file\n\n")

f.write(f"{ntot} atoms\n")
f.write(f"{nbonds} bonds\n")
f.write("0 angles\n0 dihedrals\n0 impropers\n\n")

atom_types = 1
bond_types = 1

beadType = atom_types

sigmas = [float(par['sigma_bead'])]  #list of spheres "sizes" for linked list

if par['nsolvent'] > 0:
    atom_types += 1
    solventType = atom_types
    sigmas.append(float(par['sigma_solvent']))

if par['npatch'] > 0:
    atom_types += 1
    bond_types += 2
    patchType = atom_types
    sigmas.append(float(par['sigma_patch']))

if par['ncolloids'] > 0:
    atom_types += 1
    colloidType = atom_types
    sigmas.append(float(par['sigma_colloid']))

atomList = LinkedCellList(Lbox, max(sigmas), ntot)

f.write(f"{atom_types} atom types\n{bond_types} bond types\n0 angle types\n0 dihedral types\n0 improper types\n\n")

f.write(f"0.0 {Lbox[0]} xlo xhi\n0.0 {Lbox[1]} ylo yhi\n0.0 {Lbox[2]} zlo zhi\n\n")

f.write(f"Masses\n\n{beadType} 1\n")

if par['nsolvent'] > 0:
    f.write(f"{solventType} 0.5\n")  #!Solvent mass 0.5 #TODO: read from par file

if par['npatch'] > 0:
    f.write(f"{patchType} 1\n")   #!Patch mass 1 #TODO: read from par file

if par['ncolloids'] > 0:
    f.write(f"{colloidType} 1\n")     #!Colloid mass 1 #TODO: read from par file

f.write("\n")

f.write("Atoms\n")
f.write("# atom-ID\tmol-ID\tatom-type\tx\ty\tz\n")

# -----------> POLYMERS <------------- #

coordsFirstBeads = []

minPolymerDist = par['sigma_bead']

if par['npatch'] > 0:
    minPolymerDist +=  2. * par['sigma_patch']

distBetweenBeads = 1 * par['sigma_bead'] #TODO: Define the distance between the beads depending on the chosen potential

patchyBeadsIDs = []         #List of tuples AtomID Bead + Patch to keep track of the IDs for the bonds

for p in range(0, int(par['npol'])):

    print(f'Polymer #{p+1}')
    #First bead of the polymers
    ID_start = int(1 + p * (par['ns'] + par['npatch']))

    molID = int(1 + p)

    #Random coordinate of first bead (not overlapping with other polymers)
    while True:

        coordFirstBead = [np.random.random() * par['lbox_x'],
                            np.random.random() * par['lbox_y'],
                            np.random.random() * par['lbox_z']]
        

        overlap = False

        for coord in coordsFirstBeads:

            dist_x = abs(coord[0] - coordFirstBead[0])
            dist_y = abs(coord[1] - coordFirstBead[1])

            if dist_x >= 0.5 * Lbox[0]:
                dist_x = dist_x - Lbox[0]

            if dist_y >= 0.5 * Lbox[1]:
                dist_y = dist_y - Lbox[1]

            distance2DSqrd = m.sqrt(dist_x*dist_x + dist_y*dist_y)


            if distance2DSqrd <= minPolymerDist:
                overlap = True
                break
                
        
        if overlap == False:
            break

    coordsFirstBeads.append(coordFirstBead)

    patchyBeads = []        #List of the indexes of patch beads

    for i in range(int(par['npatch'])):

        while True:

            patchyParticle = np.random.randint(par['ns'])

            if patchyParticle not in patchyBeads:

                patchyBeads.append(patchyParticle)
                break
    
    #Sort the list to have the correct molID
    patchyBeads.sort()


    #Place the beads
    for npart in range(int(par['ns'])):

        coord = coordFirstBead.copy()
        coord[2] += distBetweenBeads * npart   

        coordPBC = [c - 0.5*l for c, l in zip(coord, Lbox)]

        for ax in range(3):
            coordPBC[ax] -= Lbox[ax] * round(coordPBC[ax] / Lbox[ax])

        spheres.append(Sphere(par['sigma_bead'], coordPBC, SphereType.BEAD))
        atomList.addObjectToList(ID_start + npart - 1, coordPBC)


        f.write(f"{ID_start + npart} {molID} {beadType} {coord[0]} {coord[1]} {coord[2]}\n")

    
    #Place the patches
    patchID = int(ID_start + par['ns'])


    for nparticle in patchyBeads:
          
        patchCoord = coordFirstBead.copy()
        patchCoord[2] += distBetweenBeads * nparticle

        #Random position of the patch around the bead
        theta = np.random.rand() * 2. * m.pi

        patchCoord[0] += m.cos(theta) * 0.5 * (par['sigma_bead'] + par['sigma_patch'])
        patchCoord[1] += m.sin(theta) * 0.5 * (par['sigma_bead'] + par['sigma_patch'])

        patchCoordPBC = [c - 0.5*l for c, l in zip(patchCoord, Lbox)]

        for ax in range(3):
            patchCoordPBC[ax] -= Lbox[ax] * round(patchCoordPBC[ax] / Lbox[ax])

        spheres.append(Sphere(par['sigma_patch'], patchCoordPBC, SphereType.PATCH))

        atomList.addObjectToList(patchID - 1, [c - 0.5*l for c, l in zip(patchCoordPBC, Lbox)])

        f.write(f"{patchID} {molID} {patchType} {patchCoord[0]} {patchCoord[1]} {patchCoord[2]}\n")

        #Save the tuple AtomID of the bead + patch for the bond

        patchyBeadsIDs.append((ID_start + nparticle, patchID))

        patchID += 1


# -----------> COLLOIDS <------------- #


if par['ncolloids'] > 0:

    collID = int(1 + par['npol'] * (par['ns'] + par['npatch']))
    molID = int(1 + par['npol'])

    for nc in range(int(par['ncolloids'])):


        print(f'Colloid #{nc+1}')

        while True:

            overlap = False

            colloid = Sphere(par['sigma_colloid'], (Lbox[0] * (-0.5 + np.random.random()),
                            Lbox[1] * (-0.5 + np.random.random()),
                            Lbox[2] * (-0.5 + np.random.random())), SphereType.COLLOID)

            if atomList.overlapCheck(colloid, spheres) == False:
                break

        spheres.append(colloid)
        
        atomList.addObjectToList(collID - 1, colloid.cm)
        
        f.write(f"{collID} {molID} {colloidType} {colloid.cm[0]+0.5*Lbox[0]} {colloid.cm[1]+0.5*Lbox[1]} {colloid.cm[2]+0.5*Lbox[2]}\n")

        collID += 1
        molID += 1


# -------------> SOLVENT <--------------------------- #

if par['nsolvent'] > 0:

    solvID = int(1 + par['npol'] * (par['ns'] + par['npatch']) + par['ncolloids'])
    molID = int(1 + par['npol'] + par['ncolloids'])

    coordsSolvent = []

    for nc in range(int(par['nsolvent'])):

        print(f'Solvent #{nc+1}')

        while True:

            overlap = False

            solvent = Sphere(par['sigma_solvent'], (Lbox[0] * (-0.5 + np.random.random()),
                            Lbox[1] * (-0.5 + np.random.random()),
                            Lbox[2] * (-0.5 + np.random.random())), SphereType.SOLVENT)

            if atomList.overlapCheck(solvent, spheres) == False:
                break
        
        spheres.append(solvent)

        atomList.addObjectToList(solvID - 1, solvent.cm)
        
        f.write(f"{solvID} {molID} {solventType} {solvent.cm[0]+0.5*Lbox[0]} {solvent.cm[1]+0.5*Lbox[1]} {solvent.cm[2]+0.5*Lbox[2]}\n")

        solvID += 1
        molID += 1


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
for p in range(0, int(par['npol'])):

    #First bead of the polymer
    ID_start = int(1 + p * (par['ns'] + par['npatch']))

    for i in range(int(par['ns'] - 1)):

        atom_id = ID_start + i

        f.write(f"{ID_bond} 1 {atom_id} {atom_id + 1}\n")

        ID_bond += 1


#Bonds between the bead and the patch
for BeadID,PatchID in patchyBeadsIDs:

    f.write(f"{ID_bond} 2 {BeadID} {PatchID}\n")
    ID_bond += 1


f.close()

# -----------> .xyz for VMD visualization <------------- #

f = open("polymers.xyz", "w+")

f.write(f"{ntot}\n\n")

fakenames = {1: 'C', 2:'H', 3:'O', 4:'N'}

for sphere in spheres:

    if sphere.sphereType == SphereType.BEAD:
        letter = fakenames[beadType]
    elif sphere.sphereType == SphereType.PATCH:
        letter = fakenames[patchType]
    elif sphere.sphereType == SphereType.COLLOID:
        letter = fakenames[colloidType]
    elif sphere.sphereType == SphereType.SOLVENT:
        letter = fakenames[solventType]
    f.write(f"{letter} {sphere.cm[0]} {sphere.cm[1]} {sphere.cm[2]}\n")

f.close()