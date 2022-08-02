"""
Protein Data Bank (PDB) acquisition package (Bio.PDB.acquisition).

This package respects python module naming convention (see: https://peps.python.org/pep-0008/#package-and-module-names).

Materials:
- https://www.wwpdb.org/ftp/pdb-ftp-sites
- https://www.rcsb.org/docs/programmatic-access/file-download-services

NOTE: only export the object consumed outside of the module.
"""

# Copyright (C) 2022, Julian Maurin (julian.maurin.dev@proton.me)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.PDB.acquisition.legacy import handle_legacy_server
from Bio.PDB.acquisition.protocol import get_default_protocol
from Bio.PDB.acquisition.server import PDBServer, get_fastest_server

__all__ = [
    "get_default_protocol",
    "handle_legacy_server",
    "PDBServer",
    "get_fastest_server",
]
