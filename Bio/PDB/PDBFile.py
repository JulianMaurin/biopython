# Copyright (C) 2022, Julian Maurin (julian.maurin.dev@proton.me)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Protein Data Bank (PDB) file format module.

Provides data and logics to build PDB file path according to the file format.
"""
import dataclasses
from Bio.PDB.PDBServer import PDB_SERVERS, PDBServer, PDBServerProtocol


class BuildPDBFilePathError(Exception):
    """Error building PDB file path."""


@dataclasses.dataclass(frozen=True)
class PDBFileFormat:
    """Build PDB file path."""

    label: str
    extension: str
    directory_pattern: str
    filename_pattern: str = "{code}"
    # We assume obsolete files are available if the directory is set.
    obsolete_directory_pattern: str | None = None
    # We assume compressed files are available if the extension is set.
    compressed_extension: str | None = "gz"
    # Format may group several file for the same unit (eg. Assembly).
    # In this case the files names are numbered using the "index" variable.
    multiple_files: bool | None = None
    # Format may be available on subdomain (eg. MMTF).
    # NOTE: File format subdomain is overriding the one related to protocol.
    subdomain: str | None = None
    # Format files may be located outside of the default pdb directory.
    # In this case the directory_pattern should contains the root directory.
    # NOTE: The pdb directory is associated to the PDBServer as it may differ
    # from a server to another.
    is_located_into_pdb_directory: bool = True
    # The servers where the format is available.
    servers: tuple[PDBServer] = tuple(server for server in PDB_SERVERS.values())
    # The protocol available to download the file.
    protocols: tuple[PDBServerProtocol] = tuple(
        protocol for protocol in PDBServerProtocol
    )

    def __str__(self) -> str:
        return self.label

    def build_name(
        self, code: str, compressed: bool | None = None, index: int | None = None
    ):
        """Build file name."""
        if compressed and not self.has_compressed_file:
            raise BuildPDBFilePathError(
                f"Compressed files do not exist for this format (format: {self})."
            )

        if self.multiple_files and not index:
            raise BuildPDBFilePathError(
                f"Building file name for this format require file index (format: {self})."
            )

        # Numbered filename contains index (see. PDBFile.multiple_files) at the
        # filename or extension level.
        filename = self.filename_pattern.format(code=code, index=index)
        extension = f".{self.extension.format(index=index)}" if self.extension else ""
        extra_extension = (
            f".{self.compressed_extension}"
            if compressed and self.compressed_extension
            else ""
        )
        return f"{filename}{extension}{extra_extension}"

    def build_path(
        self,
        code: str,
        compressed: bool = True,
        obsolete: bool = False,
        index: int | None = None,
    ):
        """Build file path."""
        if obsolete and not self.has_obsolete_file:
            raise BuildPDBFilePathError(
                "Obsolete files do not exist for this format (format: {self})."
            )

        # Entries are grouped by the middle two characters of the
        # 4-character PDB identifier.
        short_code = code[1:3]

        # Structures and associated data files no longer part of
        # the archive are located in another directory.
        pattern = (
            self.obsolete_directory_pattern if obsolete else self.directory_pattern
        )

        return pattern.format(
            code=code,
            short_code=short_code,
            filename=self.build_name(code=code, compressed=compressed, index=index),
            index=index,
        )

    @property
    def has_obsolete_file(self) -> bool:
        """Obsolete files are considered to exist if the directory is set."""
        return self.obsolete_directory_pattern is not None

    @property
    def has_compressed_file(self) -> bool:
        """Compressed files are considered to exist if the extension is set."""
        return self.compressed_extension is not None


FILE_FORMATS = {
    # Legacy format for PDB entry.
    "PDB": PDBFileFormat(
        label="PDB",
        extension="ent",
        filename_pattern="pdb{code}",
        directory_pattern="/data/structures/divided/pdb/{short_code}/{filename}",
        obsolete_directory_pattern="/data/structures/obsolete/pdb/{short_code}/{filename}",
    ),
    # The biologically relevant and/or functional grouping of a particular set
    # of macromolecules (PDB format).
    "PDB_ASSEMBLY": PDBFileFormat(
        label="PDB (biological assembly)",
        extension="pdb{index}",
        directory_pattern="/data/biounit/PDB/divided/{short_code}/{filename}",
        multiple_files=True,
    ),
    # Macromolecular Crystallographic Information File.
    # Default format since 2014.
    "MMCIF": PDBFileFormat(
        label="PDBx/mmCIF",
        extension="cif",
        directory_pattern="/data/structures/divided/mmCIF/{short_code}/{filename}",
        obsolete_directory_pattern="/data/structures/obsolete/mmCIF/{short_code}/{filename}",
    ),
    # The biologically relevant and/or functional grouping of a particular set
    # of macromolecules (mmCIF format).
    "MMCIF_ASSEMBLY": PDBFileFormat(
        label="PDBx/mmCIF (biological assembly)",
        extension="cif",
        filename_pattern="{code}-assembly{index}",
        directory_pattern="/data/assemblies/mmCIF/divided/{short_code}/{filename}",
        multiple_files=True,
    ),
    # The Protein Data Bank Markup Language (PDBML) provides a
    # representation of PDB data in XML format.
    "PDBML": PDBFileFormat(
        label="PDBML/XML",
        extension="xml",
        directory_pattern="/data/structures/divided/XML/{short_code}/{filename}",
        obsolete_directory_pattern="/data/structures/obsolete/XML/{short_code}/{filename}",
    ),
    # TAR files containing a collection of best effort/minimal files in
    # the PDB format are available for some of the entries that do not
    # have legacy PDB-format files.
    "PDB_BUNDLE": PDBFileFormat(
        label="PDB bundle",
        filename_pattern="{code}-pdb-bundle",
        extension="tar",
        directory_pattern="/compatible/pdb_bundle/{short_code}/{code}/{filename}",
    ),
    # The Macromolecular Transmission Format (MMTF) is an extensible, compact,
    # self-contained binary format to efficiently transmit, load,
    # and process 3D biomolecular structural data.
    "MMTF": PDBFileFormat(
        label="MMTF",
        extension="",
        directory_pattern="/v1.0/full/{filename}",
        compressed_extension="",  # the file is compressed but there is not extension.
        servers=(PDB_SERVERS["US"],),
        subdomain="mmtf",
        protocols=(PDBServerProtocol.HTTPS,),
        is_located_into_pdb_directory=False,
    ),
}
