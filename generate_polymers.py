#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

from core.parse_parameters import parseParameters
from core.particle_system import ParticleSystem

#Read parameter file name
if len(sys.argv) != 2:
    print("Specify the name of the parameter file")
    sys.exit(1)

parFile = sys.argv[1]

partsys = ParticleSystem(parseParameters(parFile))

partsys.addPolymers()
partsys.addColloids()
partsys.addSolvent()

partsys.writeLAMMPSfile()

if partsys.par["saveXYZ"] == True:
    partsys.saveXYZfile()