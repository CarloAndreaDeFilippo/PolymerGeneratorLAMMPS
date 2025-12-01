#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
from os import listdir, getcwd
import sys

from core.parse_parameters import parseParameters
from core.make_directories import makeDirectories

def minPositionLJ(sigma: float) -> float:
    """Minimum of LJ to mimick WCA"""
    return 2**(1.0/6.0) * sigma

def maxDistanceInteractionLJ(sigma: float) -> float:
    """Maxmimum distance interaction for LJ"""
    return 2.5 * sigma

#Read parameter file name
if len(sys.argv) != 2:
    print("Specify the name of the parameter file")
    sys.exit(1)

parFile = sys.argv[1]

# Reading parameter file
par = parseParameters(parFile)

#Creation of the necessary directories
makeDirectories([getcwd() + "/configurations"])

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

#TODO: add epsilon values for all the combinations
f = open(f"{par["inputFile"]}", "w+")

f.write("######################\n#   Initialization   #\n######################\n\n")
	  
f.write("units lj\n")
f.write("atom_style bond\n")
f.write("special_bonds fene\n") #To not count the LJ potential twice for bonded beads
f.write("boundary p p p\n")
f.write("comm_modify mode single cutoff 5.0 vel yes")   #For thread communication (mpi)

f.write("\n")

f.write("#######################\n#   Atom definition   #\n#######################\n\n")


if int(par['restart']) == 1:  
    f.write("read_restart polymers.restart remap\n")
    f.write("\n")
else:
    f.write(f"read_data {par["atomFile"]} extra/special/per/atom 100\n") # da modificare per prendere parametri esterni
    f.write("\n")

f.write("################\n#   Settings   #\n################\n\n")

#*Spheres setup

#Cut and shifted LJ potential at the minimum (mimic WCA potential)
f.write(f"pair_style lj/cut {minPositionLJ(par['sigma_bead'])}\n")
f.write(f"pair_modify shift yes\n")

#LJ interaction: type1 type2 epsilon sigma cutoff

#Beads-beads interaction
f.write(f"pair_coeff {beadType} {beadType} {par['eps']} {par['sigma_bead']} {minPositionLJ(par['sigma_bead'])}\n")   #! WCA

if par['nsolvent'] > 0:
    #Bead-solvent interaction
    sigma_bead_solv = 0.5 * (par['sigma_bead'] + par['sigma_solvent'])
    f.write(f"pair_coeff {beadType} {solventType} {par['eps']} {sigma_bead_solv} {minPositionLJ(sigma_bead_solv)}\n")   #! WCA
    
    #Solvent-solvent interaction
    f.write(f"pair_coeff {solventType} {solventType} {par['eps_ss']} {par['sigma_solvent']} {maxDistanceInteractionLJ(par['sigma_solvent'])}\n") #! LJ cut

if par['npatch'] > 0:
    #Bead-patch interaction
    sigma_bead_patch = 0.5 * (par['sigma_bead'] + par['sigma_patch'])
    f.write(f"pair_coeff {beadType} {patchType} {par['eps']} {sigma_bead_patch} {minPositionLJ(sigma_bead_patch)}\n")   #! WCA
    
    #Patch-patch interaction
    f.write(f"pair_coeff {patchType} {patchType} {par['eps']} {par['sigma_patch']} {minPositionLJ(par['sigma_patch'])}\n")   #! WCA

    if par['nsolvent'] > 0:
        #Patch-solvent
        sigma_patch_solv = 0.5 * (par['sigma_patch'] + par['sigma_solvent'])
        f.write(f"pair_coeff {patchType} {solventType} {par['eps']} {sigma_patch_solv} {minPositionLJ(sigma_patch_solv)}\n")   #! WCA


if par['ncolloids'] > 0:
    #Bead-colloid interaction
    sigma_bead_coll = 0.5 * (par['sigma_bead'] + par['sigma_colloid'])
    f.write(f"pair_coeff {beadType} {colloidType} {par['eps_bc']} {sigma_bead_coll} {maxDistanceInteractionLJ(sigma_bead_coll)}\n")   #! LJ cut
    
    #Colloid-colloid interaction
    f.write(f"pair_coeff {colloidType} {colloidType} {par['eps']} {par['sigma_colloid']} {minPositionLJ(par['sigma_colloid'])}\n")    #! WCA

    if par['nsolvent'] > 0:
        #Colloid-solvent interaction
        sigma_coll_solv = 0.5 * (par['sigma_colloid'] + par['sigma_solvent'])
        f.write(f"pair_coeff {colloidType} {solventType} {par['eps_cs']} {sigma_coll_solv} {maxDistanceInteractionLJ(sigma_coll_solv)}\n")   #! LJ cut
 
    if par['npatch'] > 0:
        #Colloid-patch interaction
        sigma_coll_patch = 0.5 * (par['sigma_colloid'] + par['sigma_patch'])
        f.write(f"pair_coeff {colloidType} {patchType} {par['eps']} {sigma_coll_patch} {minPositionLJ(sigma_coll_patch)}\n")   #! WCA
    

f.write("\n")

#*Bonds setup

f.write("neigh_modify one 10000\n")

f.write("bond_style fene\n")

#FENE bond: bondID Kparameter MaxDistance LJeps sigma

#Bond between beads
f.write(f"bond_coeff 1 30.0 {1.5 * par['sigma_bead']} {par['eps']} {par['sigma_bead']}\n") 

if par['npatch'] > 0:
    #Bond between bead and patch
    sigma_bead_patch = 0.5 * (par['sigma_bead'] + par['sigma_patch'])
    f.write(f"bond_coeff 2 30.0 {1.5 * sigma_bead_patch} {par['eps']} {sigma_bead_patch}\n")

    #Bond for the crosslink between patches
    f.write(f"bond_coeff 3 30.0 {1.5 * par['sigma_patch']} {par['eps']} {par['sigma_patch']}\n")

f.write("\n")	 


#*Simulation environment setup

f.write("fix 1 all nve\n")                             #NVE integrator (combined with Langevin thermostat gives NVT)
f.write("fix 2 all langevin 1.0 1.0 1.0 16113\n")      #Langevin thermostat (Tstart Tfin dampening seed)

if par['press/berendsen'] > 0:
    f.write(f"fix 3 all press/berendsen iso {par['press_in']} {par['press_fin']} 1000.0\n") #If NPT with Berendsen barostat

#f.write("velocity all create 50.0 66254 rot yes dist gaussian\n") #Random initial velocities



f.write("#####################\n#   Calculations    #\n#####################\n")
f.write("\n")

file_number = 0

if int(par['restart']) == 1:
        file_number = int(re.search("[0-9]", [i for i in listdir("gyration") if "gyr" in i][0]).group()) + 1


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

f.write("restart 250000 polymers.restart polymers2.restart\n")  #LAMMPS restart files #TODO: read from parameter file the restart frequency
f.write("\n")
f.write(f"dump config all custom {int(par['dumpsteps'])} configurations/polymers_{file_number}.lammpstrj id mol type x y z\n")
f.write("\n")
#f.write("dump endtoend bbone_extremes xyz 500000 gyration/bbone_extremes_*.lammpstrj\n")# id type x y z\n".format(file_number))
#f.write("\n")

f.write("##########################\n#   Running simulation   #\n##########################\n")
f.write("\n")
f.write("timestep 0.001\n") #Timestep of the integrator 
f.write(f"thermo {int(par['thermosteps'])}\n")
f.write("thermo_style custom step temp etotal ke pe epair emol press vol\n")
#!thermo_style custom step temp etotal ke pe epair emol press vol c_com0[1] c_com0[2] c_com0[3]
f.write("thermo_modify flush yes\n")

#*Crosslink: bond creation during the simulation between patches
#fix fixName groupID bond/create Nevery itype jtype Rmin bondtype prob p seed iparam maxbonds newitype jparam maxbonds newjtype

#Create the crosslinks if there are patches in the system and the crosslink parameter is set to True (value != 0)
if par['npatch'] > 0 and par['crosslink'] == True:
    f.write(f"fix createBonds all bond/create 1 {patchType} {patchType} {1 * par['sigma_bead']} 3 iparam 1 2 jparam 1 2\n")

f.write("minimize 1.0e-4 1.0e-6 100 1000\n") #Little energy minimization to stabilize initial configuration

f.write(f"run {int(par['totsteps'])}\n")

f.close()
