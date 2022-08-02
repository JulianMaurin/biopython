"""Acquisition server module."""
from __future__ import annotations

import dataclasses
import functools
import time

from Bio.PDB.acquisition import exceptions
from Bio.PDB.acquisition.config import get_config
from Bio.PDB.acquisition.protocol import PDBServerProtocol
from Bio.PDB.acquisition.utils import build_socket


@dataclasses.dataclass(frozen=True)
class PDBServer:
    """PDB server immutable object."""

    code: str
    domain: str
    label: str
    pdb_directory: str

    def __str__(self) -> str:
        return f"{self.label} ({self.domain})"

    @classmethod
    def subdomain(cls, protocol: PDBServerProtocol) -> str:
        """Build subdomain URL part based on the protocol."""
        return "ftp." if protocol == PDBServerProtocol.FTP else ""

    def build_pdb_directory_url(self, protocol: PDBServerProtocol) -> str:
        """Build URL to server pdb directory."""
        return f"{protocol.url_prefix}{self.subdomain(protocol)}{self.domain}{self.pdb_directory}"


def get_pdb_servers() -> list[PDBServer]:
    """Get PDB Servers from JSON config."""
    configuration = get_config()
    try:
        return [
            PDBServer(
                code=code,
                domain=server["domain"],
                label=server["label"],
                pdb_directory=server["pdb_directory"],
            )
            for code, server in configuration["servers"].items()
        ]
    except Exception as err:
        raise exceptions.ConfigurationError(
            f"Error parsing server from configuration (configuration: {configuration})."
        ) from err


@functools.cache
def _get_server_connection_timing(
    server: PDBServer, protocol: PDBServerProtocol
) -> int | None:
    """Duration to connect to a PDB server."""
    domain = f"{server.subdomain(protocol)}{server.domain}"
    with build_socket() as server_socket:
        time_before_connect = time.time()
        try:
            server_socket.connect((domain, protocol.port))
        except Exception:
            return None
        return time.time() - time_before_connect


@functools.cache
def get_fastest_server(protocol: PDBServerProtocol) -> PDBServer:
    """Fastest server among the available servers."""
    fastest_server = None
    fastest_timing = None

    for server in get_pdb_servers():
        timing = _get_server_connection_timing(server, protocol)

        if timing and (not fastest_timing or timing < fastest_timing):
            fastest_server = server
            fastest_timing = timing

    if not fastest_server:
        raise exceptions.PDBServersConnectionError(protocol)

    return fastest_server
