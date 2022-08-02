"""Acquisition exceptions module."""
from __future__ import annotations

from Bio.PDB.PDBExceptions import PDBException


class AcquisitionException(PDBException):
    """Module top level error."""


class ConfigurationError(AcquisitionException):
    """Error related to the configuration."""
