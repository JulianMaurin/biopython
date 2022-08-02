# Copyright (C) 2022, Julian Maurin (julian.maurin.dev@proton.me)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Protein Data Bank (PDB) server module.

Materials:
- https://www.wwpdb.org/ftp/pdb-ftp-sites
- https://www.rcsb.org/docs/programmatic-access/file-download-services
"""

from __future__ import annotations

import dataclasses
import enum
import functools
import re
import time

from Bio.PDB.acquisition.utils import build_socket

SERVER_REGEX = re.compile(
    r"^(?P<protocol>(http|ftp)s?):\/{2}(ftp\.)?(?P<domain>.*\..*?)(\/)?$"
)


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


PDB_SERVERS = [
    PDBServer(
        code="WW", domain="wwpdb.org", label="Worldwide", pdb_directory="/pub/pdb"
    ),
    PDBServer(
        code="US", domain="rcsb.org", label="United States", pdb_directory="/pub/pdb"
    ),
    PDBServer(
        code="UK",
        domain="ebi.ac.uk",
        label="United Kingdom",
        pdb_directory="/pub/databases/pdb",
    ),
    PDBServer(code="JP", domain="pdbj.org", label="Japan", pdb_directory="/pub/pdb"),
]


class PDBServerProtocol(enum.IntEnum):
    """Available protocol to access to PDB servers."""

    HTTPS = 443
    FTP = 21

    @property
    def port(self) -> int:
        """Protocol port."""
        return self.value

    @property
    def url_prefix(self) -> str:
        """Build protocol url prefix."""
        return f"{self.name.lower()}://"

    def __str__(self) -> str:
        return f"{self.name} ({self.value})"


class PDBServersConnectionError(Exception):
    """Unable to connect to any PDB server."""

    def __init__(self, protocol: PDBServerProtocol) -> None:
        """Build PDBServersConnectionError message."""
        super().__init__(
            f"Unable to connect to any PDB servers (protocol: {protocol.name}, servers: {', '.join([str(server) for server in PDB_SERVERS])})."
        )


class UnsupportedServerError(Exception):
    """PDB Server is not supported."""

    def __init__(self, server: str) -> None:
        """Build UnsupportedServerError message."""
        super().__init__(
            f"PDB server is not supported (server: {server}, supported servers: {', '.join([str(server) for server in PDB_SERVERS])})."
        )


class UnsupportedProtocolError(Exception):
    """Protocol is not supported."""

    def __init__(self, protocol: str) -> None:
        """Build UnsupportedProtocolError message."""
        super().__init__(
            f"Protocol is not supported (protocol: {protocol}, supported protocol: {', '.join([str(_protocol) for _protocol in PDBServerProtocol])})."
        )


@functools.cache
def get_server_connection_timing(
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

    for server in PDB_SERVERS:
        timing = get_server_connection_timing(server, protocol)

        if timing and (not fastest_timing or timing < fastest_timing):
            fastest_server = server
            fastest_timing = timing

    if not fastest_server:
        raise PDBServersConnectionError(protocol)

    return fastest_server


def handle_legacy_server(server: str) -> tuple(PDBServer, PDBServerProtocol):
    """Provide support for legacy server declaration (as string)."""
    try:
        match = SERVER_REGEX.match(server.lower())
        protocol_str = match["protocol"]
        domain_str = match["domain"]
    except Exception:
        raise UnsupportedServerError(server=server)

    protocol_obj = None
    for available_protocol in PDBServerProtocol:
        if available_protocol.name.lower() == protocol_str:
            protocol_obj = available_protocol
            break
    if not protocol_obj:
        raise UnsupportedProtocolError(protocol=protocol_str)

    server_obj = None
    for available_server in PDB_SERVERS:
        if available_server.domain == domain_str:
            server_obj = available_server
            break
    if not server_obj or not protocol_obj:
        raise UnsupportedServerError(server=server)

    return (server_obj, protocol_obj)
