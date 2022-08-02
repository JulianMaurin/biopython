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
import time

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
