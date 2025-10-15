# PolymerGeneratorLAMMPS

A Python-based tool to automatically generate [LAMMPS](https://github.com/lammps/lammps) configurations and input files for polymer-solvent systems with optional crosslinking and crowding effects.

<img src="/docs/systemExample.png" width=50% height=50%>

## How to use the script

1. Edit simulation parameters in `par.txt` (e.g., number of polymers, solvent, interaction strengths).
2. Generate the configuration file:

```
python3 generate_configuration.py
```

3. Create the corresponding LAMMPS input file:
   
```
python3 generate_input.py
```
4. Run the simulation:

```
lmp_serial -in in.polymers
```

or, in parallel:

```
mpirun -np NUM_THREADS lmp_mpi -in in.polymers
```

## Simulation details

### Configuration file
The `generate_polymers.py` script creates a LAMMPS configuration file (`polymers.dat`) containing:
* `npol`: number of self-avoiding polymers of equal length `ns` (WCA interaction between the beads)
* `npatch`: number of additional patchy spheres on each polymer, placed randomly along the chain; these can be used to simulate crosslinking
* `nsolvent`: number of solvent spheres placed randomly in the box. These interact with polymer beads via WCA and with other solvent atoms via Lennard-Jones potentials.
* `ncolloids`: number of colloidal spheres representing additional crowders (LJ interaction with both polymer beads and solvent)

The sizes of all spheres (polymer beads, solvent, colloids, etc.) can be tuned using the various `sigma_` parameters in the `par.txt` file.

A cell linked-list is implemented to accelerate the random placement of spheres within the simulation box.

### Input file
The `generate_input.py` script generates the input file for LAMMPS:
* The simulation can be either NVT (Langevin thermostat for Brownian dynamics) or NPT (Berendsen barostat if `press/berendsen` not 0)
* `totsteps`: total number of MD steps of the simulation
* `thermosteps`: frequency of thermo output to the terminal
* `dumpsteps`: frequency of saving trajectory snapshots to the ```.lammpstrj``` file
* `crosslink`: switch on/off of the crosslinker interaction. If non-zero, two patches that approach within a threshold distance will form a permanent FENE bond.


## Output files

The scripts will give three files:
* `generate_polymers.py` &#8594; `polymers.dat` (LAMMPS configuration) and `polymers.xyz` (viewable with [VMD](https://www.ks.uiuc.edu/Research/vmd/))
* `generate_input.py` &#8594; `in.polymers` (LAMMPS input script)
