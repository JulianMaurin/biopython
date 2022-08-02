"""Acquisition exceptions module."""
from __future__ import annotations

from typing import TYPE_CHECKING

from Bio.PDB.acquisition.server import get_pdb_servers
from Bio.PDB.PDBExceptions import PDBException

if TYPE_CHECKING:
    from Bio.PDB.acquisition.protocol import PDBServerProtocol


class AcquisitionException(PDBException):
    """Module top level error."""


class ConfigurationError(AcquisitionException):
    """Error related to the configuration."""


class PDBServersConnectionError(AcquisitionException):
    """Unable to connect to any PDB server."""

    def __init__(self, protocol: PDBServerProtocol) -> None:
        """Build PDBServersConnectionError message."""
        super().__init__(
            f"Unable to connect to any PDB servers (protocol: {protocol.name}, servers: {', '.join([str(server) for server in get_pdb_servers()])})."
        )


class UnsupportedServerError(AcquisitionException):
    """PDB Server is not supported."""

    def __init__(self, server: str) -> None:
        """Build UnsupportedServerError message."""
        super().__init__(
            f"PDB server is not supported (server: {server}, supported servers: {', '.join([str(server) for server in get_pdb_servers()])})."
        )


class UnsupportedProtocolError(AcquisitionException):
    """Protocol is not supported."""

    def __init__(self, protocol: str) -> None:
        """Build UnsupportedProtocolError message."""
        from Bio.PDB.acquisition.protocol import PDBServerProtocol

        super().__init__(
            f"Protocol is not supported (protocol: {protocol}, supported protocol: {', '.join([str(_protocol) for _protocol in PDBServerProtocol])})."
        )
