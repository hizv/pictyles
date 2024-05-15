#!/usr/bin/env python3
import sys
from os import path

sys.path.append(path.dirname(path.dirname(__file__)))

from charmpic.pic import CCSInterface, Species, Simulation


interface = CCSInterface(sys.argv[1], 1234)

electron = Species(
    name="electron",
    init_count=16000,
    mass=1,
    charge=10,
    position_init="linear",
    pos_dist_params=[1, 2],
)

sim = Simulation(
    interface=interface,
    odf=16,
    sim_box_shape=[5, 5],
    iterations=4000,
    boundary_conditions="periodic",
    species=[electron],
    geometry="cartesian2D",
    lb_freq=1200,
    migrate_freq=1,
)

sim.run()
