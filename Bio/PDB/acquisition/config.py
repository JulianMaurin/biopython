"""
Acquisition configuration module, loading the JSON configuration file.

NOTE:   Do not import any biopython at the module
        level to avoid circular import.
"""
from __future__ import annotations

import functools
import json
import os
import typing


CONFIG_FILE_PATH = os.path.join(os.path.dirname(__file__), "config.json")


class Default(typing.TypedDict):
    """Configuration type (default)."""

    protocol: str


class Server(typing.TypedDict):
    """Configuration type (server)."""

    domain: str
    label: str
    pdb_directory: str


class AcquisitionConfig(typing.TypedDict):
    """Configuration type (top level)."""

    default: Default
    servers: dict[str, Server]


@functools.cache
def get_config() -> AcquisitionConfig:
    """Get configuration from JSON file."""
    from Bio.PDB.acquisition.exceptions import ConfigurationError

    try:
        with open(CONFIG_FILE_PATH, "r") as file_steam:
            return json.load(file_steam)
    except Exception as err:
        raise ConfigurationError(
            f"Error loading configuration (path: {CONFIG_FILE_PATH})."
        ) from err
