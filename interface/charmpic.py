#!/usr/bin/env python3
import struct
import sys
from pyccs import Server


def to_bytes(value, dtype="I"):
    return struct.pack(dtype, value)


def from_bytes(bvalue, dtype="I"):
    return struct.unpack(dtype, bvalue)[0]


DISTRIBUTION = {"linear": 0, "sine": 1, "geometric": 2}
BOUNDARY = {"periodic": 0, "stop": 1}


class Handlers(object):
    connection_handler = b"pic_connect"
    disconnection_handler = b"pic_disconnect"
    operation_handler = b"pic_operation"
    fetch_handler = b"pic_fetch"
    delete_handler = b"pic_delete"
    exit_handler = b"pic_exit"
    sync_handler = b"pic_sync"
    create_handler = b"pic_create"


class CCSInterface(object):
    def __init__(self, server_ip, server_port):
        self.server = Server(server_ip, server_port)
        self.server.connect()
        self.client_id = self.send_command(Handlers.connection_handler, "")

    def __del__(self):
        # self.disconnect()
        pass

    def disconnect(self):
        cmd = to_bytes(self.client_id, "B")
        self.send_command_async(Handlers.disconnection_handler, cmd)

    def send_command_raw(self, handler, msg, reply_size):
        self.server.send_request(handler, 0, msg)
        return self.server.receive_response(reply_size)

    def send_command(self, handler, msg, reply_size=1, reply_type="B"):
        return from_bytes(self.send_command_raw(handler, msg, reply_size), reply_type)

    def send_command_async(self, handler, msg):
        self.server.send_request(handler, 0, msg)


class Species(object):
    def __init__(
        self,
        name: str,
        particles_per_cell: int,
        mass: float = 1.0,
        charge: float = 1.0,
        position_init: str = "linear",
    ):
        self.name = name
        self.particles_per_cell = particles_per_cell
        self.position_init = position_init  # [linear, sine, geometric]
        self.mass = mass
        self.charge = charge


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

    def run(self):
        """
        grid_shape
        odf
        sim_box_shape
        time_delta
        """
        # cmd = to_bytes(len(self.grid_shape), "B")
        # for dim in self.grid_shape:
        #     cmd += to_bytes(dim, "I")

        cmd = to_bytes(len(self.sim_box_shape), "B")
        for dim in self.sim_box_shape:
            cmd += to_bytes(dim, "I")

        cmd += to_bytes(self.odf, "B")
        cmd += to_bytes(self.time_delta, "f")
        cmd += to_bytes(self.iterations, "I")
        cmd += to_bytes(BOUNDARY[self.boundary_conditions], "B")

        cmd += to_bytes(self.species[0].particles_per_cell, "I")
        cmd += to_bytes(DISTRIBUTION[self.species[0].position_init], "B")

        self.interface.send_command(Handlers.create_handler, cmd)


if __name__ == "__main__":
    interface = CCSInterface(sys.argv[1], 1234)

    electron = Species(name="electron", particles_per_cell=10000, mass=0.5, charge=0.1)

    sim = Simulation(
        interface=interface,
        odf=16,
        sim_box_shape=[100, 100],
        iterations=1000,
        boundary_conditions="periodic",
        species=[electron],
    )

    sim.run()
