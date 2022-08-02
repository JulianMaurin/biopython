"""Acquisition legacy module.

Provides functions to ensure backward compatibility.
"""
from __future__ import annotations

import re

from Bio.PDB.acquisition import exceptions
from Bio.PDB.acquisition.protocol import PDBServerProtocol
from Bio.PDB.acquisition.server import PDBServer, get_pdb_servers

_SERVER_REGEX = re.compile(
    r"^(?P<protocol>(http|ftp)s?):\/{2}(ftp\.)?(?P<domain>.*\..*?)(\/)?$"
)


def handle_legacy_server(server: str) -> tuple(PDBServer, PDBServerProtocol):
    """Provide support for legacy server declaration (as string)."""
    try:
        match = _SERVER_REGEX.match(server.lower())
        protocol_str = match["protocol"]
        domain_str = match["domain"]
    except Exception:
        raise exceptions.UnsupportedServerError(server=server)

    protocol_obj = None
    for available_protocol in PDBServerProtocol:
        if available_protocol.name.lower() == protocol_str:
            protocol_obj = available_protocol
            break
    if not protocol_obj:
        raise exceptions.UnsupportedProtocolError(protocol=protocol_str)

    server_obj = None
    for available_server in get_pdb_servers():
        if available_server.domain == domain_str:
            server_obj = available_server
            break
    if not server_obj or not protocol_obj:
        raise exceptions.UnsupportedServerError(server=server)

    return (server_obj, protocol_obj)
