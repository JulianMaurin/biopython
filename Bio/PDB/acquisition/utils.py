"""Acquisition utilities module."""

import contextlib
import socket


@contextlib.contextmanager
def build_socket():
    """Context manager auto closing socket."""
    server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    server_socket.settimeout(1)
    try:
        yield server_socket
    finally:
        server_socket.close()
