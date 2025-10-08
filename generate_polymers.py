#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import errno
import math as m
from os import makedirs, path, getcwd
from enum import Enum, auto


#Directories che ci servono

dir_list = ["configurations", "gyration",  "com"]

main_dir = getcwd() + '/'

for i in dir_list:

    try:
        makedirs(main_dir + i, 0o777)
    except OSError as exc:
        if exc.errno == errno.EEXIST and path.isdir(main_dir + i):
            pass
        else: raise


class SphereType(Enum):
    GENERIC = auto()
    BEAD = auto()
    PATCH = auto()
    COLLOID = auto()
    SOLVENT = auto()


class Sphere:

    def __init__(self, sigma=1, cm=(0., 0., 0.), sphereType=SphereType.GENERIC):
        
        self.sigma = sigma
        self.cm = cm
        self.sphereType = sphereType

    def distance(self, other, Lbox):

        dist = [0., 0., 0.]

        for ax in range(3):
            dist[ax] = self.cm[ax] - other.cm[ax]
            dist[ax] = dist[ax] - Lbox[ax] * round(dist[ax] / Lbox[ax])
        
        return m.sqrt(dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2])

    def overlap(self, other, Lbox):

        if self.distance(other, Lbox) <= 0.5 * (self.sigma + other.sigma):

            return True

        return False



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

if par['nsolvent'] > 0:
    atom_types += 1
    solventType = atom_types

if par['npatch'] > 0:
    atom_types += 1
    bond_types += 2
    patchType = atom_types

if par['ncolloids'] > 0:
    atom_types += 1
    colloidType = atom_types

f.write(f"{atom_types} atom types\n{bond_types} bond types\n0 angle types\n0 dihedral types\n0 improper types\n\n")

f.write(f"0.0 {Lbox[0]} xlo xhi\n0.0 {Lbox[1]} ylo yhi\n0.0 {Lbox[2]} zlo zhi\n\n")

f.write(f"Masses\n\n{beadType} 1\n")

if par['nsolvent'] > 0:
    f.write(f"{solventType} 0.5\n")  #!Massa solvente 0.5

if par['npatch'] > 0:
    f.write(f"{patchType} 1\n")   #!Massa patch 1

if par['ncolloids'] > 0:
    f.write(f"{colloidType} 1\n")     #! Massa colloide 1

f.write("\n")

f.write("Atoms\n")
f.write("# atom-ID\tmol-ID\tatom-type\tx\ty\tz\n")

# Generazione posizioni polimeri

coordsFirstBeads = []


minPolymerDist = 3 * par['sigma']
distBetweenBeads = 1 * par['sigma']

patchyBeadsIDs = []         #lista tuple AtomID Bead + Patch per definire i bond

for p in range(0, int(par['npol'])):

    print(f'Polymer #{p+1}')
    #primo atomo del polimero
    ID_start = int(1 + p * (par['ns'] + par['npatch']))

    molID = int(1 + p)

    #Coord primo bead + distanza minima fra polimeri
    while True:

        coordFirstBead = [np.random.random() * par['lbox_x'],
                            np.random.random() * par['lbox_y'],
                            np.random.random() * par['lbox_z']]
        

        overlap = False

        for coord in coordsFirstBeads:

            dist_x = abs(coord[0] - coordFirstBead[0])
            dist_y = abs(coord[1] - coordFirstBead[1])

            if dist_x >= 0.5 * par['lbox_x']:
                dist_x = dist_x - par['lbox_x']

            if dist_y >= 0.5 * par['lbox_y']:
                dist_y = dist_y - par['lbox_y']

            distance2DSqrd = m.sqrt(dist_x*dist_x + dist_y*dist_y)


            if distance2DSqrd <= minPolymerDist:
                overlap = True
                break
                
        
        if overlap == False:
            break

    coordsFirstBeads.append(coordFirstBead)

    patchyBeads = []        #Lista degli indici dei patchy beads

    for i in range(int(par['npatch'])):

        while True:

            patchyParticle = np.random.randint(par['ns'])

            if patchyParticle not in patchyBeads:

                patchyBeads.append(patchyParticle)
                break
    
    #Sort della lista per avere i molID corretti
    patchyBeads.sort()


    #Metto sfere dure
    for npart in range(int(par['ns'])):

        coord = coordFirstBead.copy()
        coord[2] += distBetweenBeads * npart   #! Definire la distanza fra beads iniziale in funzione del potenziale
        
        spheres.append(Sphere(par['sigma'], coord, SphereType.BEAD))

        f.write(f"{ID_start + npart} {molID} {beadType} {coord[0]} {coord[1]} {coord[2]}\n")

    
    #Metto patch
    patchID = int(ID_start + par['ns'])


    for nparticle in patchyBeads:
          
        patchCoord = coordFirstBead.copy()
        patchCoord[2] += distBetweenBeads * nparticle

        theta = np.random.rand() * 2. * m.pi

        patchCoord[0] += m.cos(theta) * 0.5 * (par['sigma'] + par['sigma_patch'])        #! Da scegliere dove piazzare la patch
        patchCoord[1] += m.sin(theta) * 0.5 * (par['sigma'] + par['sigma_patch'])

        spheres.append(Sphere(par['sigma_patch'], patchCoord, SphereType.PATCH))

        f.write(f"{patchID} {molID} {patchType} {patchCoord[0]} {patchCoord[1]} {patchCoord[2]}\n")

        #Salvo la tupla degli AtomID del Bead e della patch per il bond

        patchyBeadsIDs.append((ID_start + nparticle, patchID))

        patchID += 1




# -----------> COLLOIDI <------------- #


if par['ncolloids'] > 0:

    collID = int(1 + par['npol'] * (par['ns'] + par['npatch']))
    molID = int(1 + par['npol'])

    for nc in range(int(par['ncolloids'])):


        print(f'Colloid #{nc+1}')

        while True:

            overlap = False

            colloid = Sphere(par['sigma_colloids'], (np.random.random() * Lbox[0],
                            np.random.random() * Lbox[1],
                            np.random.random() * Lbox[2]), SphereType.COLLOID)

            for sphere in spheres:

                if sphere.overlap(colloid, Lbox):
                    overlap = True
                    break

            if overlap == False:
                break

        
        f.write(f"{collID} {molID} {colloidType} {colloid.cm[0]} {colloid.cm[1]} {colloid.cm[2]}\n")

        collID += 1
        molID += 1

        spheres.append(colloid)


# -------------> SOLVENTE <--------------------------- #

if par['nsolvent'] > 0:

    solvID = int(1 + par['npol'] * (par['ns'] + par['npatch']) + par['ncolloids'])
    molID = int(1 + par['npol'] + par['ncolloids'])

    coordsSolvent = []

    for nc in range(int(par['nsolvent'])):

        print(f'Solvent #{nc+1}')

        while True:

            overlap = False

            solvent = Sphere(par['sigma_solvent'], (np.random.random() * Lbox[0],
                            np.random.random() * Lbox[1],
                            np.random.random() * Lbox[2]), SphereType.SOLVENT)

            for sphere in spheres:

                if sphere.overlap(solvent, Lbox):
                    overlap = True
                    break

            if overlap == False:
                break

        
        f.write(f"{solvID} {molID} {solventType} {solvent.cm[0]} {solvent.cm[1]} {solvent.cm[2]}\n")

        solvID += 1
        molID += 1

        spheres.append(solvent)


f.write("\n")


#######VELOCITÀ

#f.write("Velocities\n\n")

#[f.write("{0} {1} {2} {3}\n".format(i, np.random.random() * 2., np.random.random()  * 2., np.random.random()  * 2. - 1.)) for i in range(1, ntot + 1)]

#f.write("\n")

# Generazione bond IDs
f.write("Bonds\n#bond-ID\tbond-type\tfirst-atom\tsecond-atom\n")

ID_bond = 1

#Bonds fra Beads
for p in range(0, int(par['npol'])):

    #primo bead del polimero
    ID_start = int(1 + p * (par['ns'] + par['npatch']))

    for i in range(int(par['ns'] - 1)):

        atom_id = ID_start + i

        f.write(f"{ID_bond} 1 {atom_id} {atom_id + 1}\n")

        ID_bond += 1


#Bonds fra Bead e patch
for BeadID,PatchID in patchyBeadsIDs:

    f.write(f"{ID_bond} 2 {BeadID} {PatchID}\n")
    ID_bond += 1


f.close()

#! --------------> file xyz per vmd <----------------

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