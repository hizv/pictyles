#!/usr/bin/env python3
from charmpic.interface import CCSInterface, Handlers, to_bytes

DISTRIBUTION = {"linear": 0, "sine": 1, "geometric": 2}
BOUNDARY = {"periodic": 0, "reflective": 1}
GEOMETRY = {"cartesian2D": 0, "cartesian3D": 1}


class Species(object):
    def __init__(
        self,
        name: str,
        init_count: int,
        mass: float = 1.0,
        charge: float = 1.0,
        position_init: str = "linear",
        pos_dist_params: list[float] = [1, 1],
    ):
        self.name = name
        self.init_count = init_count
        self.position_init = position_init  # [linear, sine, geometric]
        self.mass = mass
        self.charge = charge
        self.pos_dist_params = pos_dist_params


class Simulation(object):
    def __init__(
        self,
        interface: CCSInterface,
        # grid_shape: list[int],
        odf: int,
        sim_box_shape: list[int],
        species: list[Species],
        time_delta: float = 0.014,
        sim_time=None,
        iterations: int = 100,
        boundary_conditions: str = "periodic",
        geometry: str = "cartesian2D",
    ):
        # self.grid_shape = grid_shape
        self.odf = odf
        self.time_delta = time_delta
        if sim_time:
            self.sim_time = sim_time
            self.iterations = sim_time // time_delta
        else:
            self.iterations = iterations
            self.sim_time = iterations * time_delta
        self.sim_box_shape = sim_box_shape
        self.boundary_conditions = boundary_conditions
        self.interface = interface
        self.species = species
        self.geometry = geometry

    def run(self):
        """
        sim_box_shape
        odf
        time_delta (dt)
        iterations
        boundary_conditions

        #foreach species:
            init_count
            initial position_distribution
            mass
            charge
        """
        cmd = to_bytes(len(self.sim_box_shape), "B")
        for dim in self.sim_box_shape:
            cmd += to_bytes(dim, "I")

        cmd += to_bytes(self.odf, "B")
        cmd += to_bytes(self.time_delta, "f")
        cmd += to_bytes(self.iterations, "I")
        cmd += to_bytes(BOUNDARY[self.boundary_conditions], "B")
        cmd += to_bytes(GEOMETRY[self.geometry], "B")

        cmd += to_bytes(self.species[0].init_count, "I")
        cmd += to_bytes(DISTRIBUTION[self.species[0].position_init], "B")
        for param in self.species[0].pos_dist_params:
            cmd += to_bytes(param, "f")
        for _ in range(2 - len(self.species[0].pos_dist_params)):
            cmd += to_bytes(0.0, "f")
        cmd += to_bytes(self.species[0].mass, "f")
        cmd += to_bytes(self.species[0].charge, "f")

        self.interface.send_command(Handlers.create_handler, cmd)
