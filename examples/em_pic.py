#!/usr/bin/env python3
import sys
from os import path
sys.path.append(path.dirname(path.dirname(__file__)))

from charmpic.pic import CCSInterface, Species, Simulation


interface = CCSInterface(sys.argv[1], 1234)

electron = Species(name="electron", init_count=100000,
                   mass=0.5, charge=0.1, position_init="linear", pos_dist_params=[0.1, 0.5])

sim = Simulation(
    interface=interface,
    odf=16,
    sim_box_shape=[100, 100],
    iterations=1000,
    boundary_conditions="periodic",
    species=[electron],
)

sim.run()
