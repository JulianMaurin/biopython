"""Acquisition protocols module.

Materials:
- https://www.cloudflare.com/learning/network-layer/what-is-a-protocol/
"""
from __future__ import annotations

import enum

from Bio.PDB.acquisition.config import get_config
from Bio.PDB.acquisition import exceptions


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


def get_default_protocol() -> PDBServerProtocol:
    """Get default protocol from config."""
    configuration = get_config()

    try:
        protocol = configuration["default"]["protocol"]
    except Exception as err:
        raise exceptions.ConfigurationError(
            f"Error parsing default protocol from configuration (configuration: {configuration})."
        ) from err

    try:
        return PDBServerProtocol[protocol]
    except KeyError:
        raise exceptions.UnsupportedProtocolError(protocol)
