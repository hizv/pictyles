#!/usr/bin/env python3
import struct
from pyccs import Server


def to_bytes(value, dtype="I"):
    return struct.pack(dtype, value)


def from_bytes(bvalue, dtype="I"):
    return struct.unpack(dtype, bvalue)[0]


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
