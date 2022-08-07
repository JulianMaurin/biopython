# Copyright (C) 2022, Julian Maurin (julian.maurin.dev@proton.me)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Protein Data Bank (PDB) server module.

Provides logics and data to interact with PDB Server.

Materials:
- https://www.wwpdb.org/ftp/pdb-ftp-sites
- https://www.rcsb.org/docs/general-help/glossary
- https://www.rcsb.org/docs/programmatic-access/file-download-services
"""

from __future__ import annotations

import contextlib
import dataclasses
import enum
import functools
import re
import socket
import time
from pathlib import Path
from typing import TYPE_CHECKING
from urllib import request

if TYPE_CHECKING:
    from Bio.PDB.PDBFile import PDBFileFormat


def install_urllib_opener() -> None:
    """HTTP request to PDB server without user agent header may be rejected."""
    opener = request.build_opener()
    opener.addheaders = [("User-agent", "biopython")]
    request.install_opener(opener)


install_urllib_opener()


class PDBServerProtocol(enum.IntEnum):
    """Available protocol to access to PDB servers."""

    HTTPS = 443
    FTP = 21

    @property
    def port(self) -> int:
        """Protocol port."""
        return self.value

    @property
    def method(self) -> str:
        """Protocol method."""
        return self.name.lower()

    def __str__(self) -> str:
        return f"{self.name} ({self.value})"


@dataclasses.dataclass(frozen=True)
class PDBServer:
    """PDB server immutable object."""

    domain: str
    label: str
    pdb_directory: str = "/pub/pdb"
    ftp_subdomain: str = "ftp"
    http_subdomain: str = "files"

    def __str__(self) -> str:
        return f"{self.label} ({self.domain})"

    def subdomain(self, protocol: PDBServerProtocol, file_format: PDBFileFormat) -> str:
        """Build URL subdomain based on the protocol or on the file format."""
        if protocol not in file_format.protocols:
            raise UnsupportedProtocolError(
                f"File format is not available using this protocol  (format: {file_format}, protocol: {protocol}, server: {self})."
            )
        # The subdomain from the format is overriding the one from the protocol.
        if file_format.subdomain:
            return file_format.subdomain
        if protocol == PDBServerProtocol.FTP:
            return self.ftp_subdomain
        if protocol == PDBServerProtocol.HTTPS:
            return self.http_subdomain
        raise UnsupportedProtocolError(
            f"Unable to determine the subdomain corresponding to this protocol (protocol: {protocol}, server: {self})."
        )

    def build_full_domain(
        self, protocol: PDBServerProtocol, file_format: PDBFileFormat
    ) -> str:
        """Aggregate domain and subdomain."""
        return f"{self.subdomain(protocol, file_format)}.{self.domain}"

    def build_full_domain_legacy(self, protocol: PDBServerProtocol) -> str:
        """Build full domain without specifying the file format."""
        if protocol == PDBServerProtocol.FTP:
            return f"{self.ftp_subdomain}.{self.domain}"
        elif protocol == PDBServerProtocol.HTTPS:
            return f"{self.http_subdomain}.{self.domain}"
        else:
            raise UnsupportedProtocolError(
                f"Unable to determine the subdomain corresponding to this protocol (protocol: {protocol}, server: {self})."
            )

    def build_pdb_directory_url(self, protocol: PDBServerProtocol) -> str:
        """Build URL to server pdb directory (DEPRECATED).

        This is not taking account the file format.
        As expected by the legacy call.
        """
        return f"{protocol.method}://{self.build_full_domain_legacy(protocol)}{self.pdb_directory}"

    def build_file_url(
        self,
        code: str,
        file_format: PDBFileFormat,
        protocol: PDBServerProtocol,
        compressed: bool = False,
        obsolete: bool = False,
        index: int | None = None,
    ):
        """Build URL to the file on the server."""
        return (
            f"{protocol.method}://"
            f"{self.build_full_domain(protocol, file_format)}"
            f"{self.pdb_directory if file_format.is_located_into_pdb_directory else ''}"
            f"{file_format.build_path(code=code, compressed=compressed, obsolete=obsolete, index=index)}"
        )

    def retrieve_file(
        self,
        path: Path,
        code: str,
        file_format: PDBFileFormat,
        protocol: PDBServerProtocol,
        compressed: bool = True,
        obsolete: bool = False,
        index: int | None = None,
    ):
        """Download file from the server and save it locally."""
        url = self.build_file_url(
            code=code,
            file_format=file_format,
            protocol=protocol,
            compressed=compressed,
            obsolete=obsolete,
            index=index,
        )
        request.urlretrieve(url, str(path))

    def __del__(self):
        request.urlcleanup()


DEFAULT_PROTOCOL = PDBServerProtocol.FTP

PDB_SERVERS = {
    "WW": PDBServer(domain="wwpdb.org", label="Worldwide"),
    "US": PDBServer(domain="rcsb.org", label="United States"),
    "UK": PDBServer(
        domain="ebi.ac.uk",
        label="United Kingdom",
        pdb_directory="/pub/databases/pdb",
        http_subdomain="ftp",
    ),
    "JP": PDBServer(domain="pdbj.org", label="Japan", http_subdomain="data.pdbjbk1"),
}


class PDBServersConnectionError(Exception):
    """Unable to connect to any PDB server."""

    def __init__(self, protocol: PDBServerProtocol) -> None:
        """Build PDBServersConnectionError message."""
        servers = ", ".join([str(server) for server in PDB_SERVERS.values()])
        super().__init__(
            f"Unable to connect to any PDB servers (protocol: {protocol.name}, servers: {servers})."
        )


class UnsupportedServerError(Exception):
    """PDB Server is not supported."""

    def __init__(self, server: str) -> None:
        """Build UnsupportedServerError message."""
        servers = ", ".join([str(server) for server in PDB_SERVERS.values()])
        super().__init__(
            f"PDB server is not supported (server: {server}, supported servers: {servers})."
        )


class UnsupportedProtocolError(Exception):
    """Protocol is not supported."""


@functools.cache
def get_server_connection_timing(
    server: PDBServer, protocol: PDBServerProtocol
) -> int | None:
    """Duration to connect to a PDB server."""
    with _build_socket() as server_socket:
        time_before_connect = time.time()
        try:
            server_socket.connect(
                (server.build_full_domain_legacy(protocol), protocol.port)
            )
        except Exception:
            return None
        return time.time() - time_before_connect


@contextlib.contextmanager
def _build_socket():
    """Context manager auto closing socket."""
    server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    server_socket.settimeout(1)
    try:
        yield server_socket
    finally:
        server_socket.close()


@functools.cache
def get_fastest_server(protocol: PDBServerProtocol) -> PDBServer:
    """Fastest server among the available servers."""
    fastest_server = None
    fastest_timing = None

    for server in PDB_SERVERS.values():
        timing = get_server_connection_timing(server, protocol)

        if timing and (not fastest_timing or timing < fastest_timing):
            fastest_server = server
            fastest_timing = timing

    if not fastest_server:
        raise PDBServersConnectionError(protocol)

    return fastest_server


SERVER_REGEX = re.compile(
    r"^(?P<protocol>\w*):\/\/((?P<subdomain>.*?)\.)?(?P<domain>wwpdb\.org|rcsb\.org|ebi\.ac\.uk|pdbj\.org)(\/)?$"
)


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
        supported_protocols = ", ".join(
            [str(_protocol) for _protocol in PDBServerProtocol]
        )
        raise UnsupportedProtocolError(
            (
                "Protocol is not supported "
                f"(protocol: {protocol_str}, supported protocols: {supported_protocols})."
            )
        )

    server_obj = None
    for available_server in PDB_SERVERS.values():
        if available_server.domain == domain_str:
            server_obj = available_server
            break
    if not server_obj or not protocol_obj:
        raise UnsupportedServerError(server=server)

    return (server_obj, protocol_obj)
