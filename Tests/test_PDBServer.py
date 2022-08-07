# Copyright (C) 2022, Julian Maurin (julian.maurin.dev@proton.me)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Testing PDB server utils."""

import urllib.error
import os
import pathlib
import tempfile
import unittest
import unittest.mock
import uuid
from Bio.PDB.PDBFile import FILE_FORMATS, PDBFileFormat

import requires_internet

# We want to test this module:
from Bio.PDB import PDBServer

requires_internet.check()


class TestGetFastestDefaultServer(unittest.TestCase):
    """Test PDBServer.get_fastest_server."""

    def setUp(self):
        super().setUp()
        PDBServer.get_fastest_server.cache_clear()
        PDBServer.get_server_connection_timing.cache_clear()

    @unittest.mock.patch("Bio.PDB.PDBServer.time.time", side_effect=[0, 5, 5, 6])
    @unittest.mock.patch(
        "Bio.PDB.PDBServer.PDB_SERVERS",
        {
            "WW": PDBServer.PDB_SERVERS["WW"],
            "US": PDBServer.PDB_SERVERS["US"],
        },
    )
    @unittest.mock.patch("Bio.PDB.PDBServer.socket.socket.connect", return_value=None)
    def test_timing(self, _, __):
        """Check that the returned server is the fastest.

        Mocking the WW server duration to 5s (5-0) and the US to 1s (6-5).
        """
        self.assertEqual(
            PDBServer.get_fastest_server(PDBServer.PDBServerProtocol.FTP),
            PDBServer.PDB_SERVERS["US"],
        )

    @unittest.mock.patch(
        "Bio.PDB.PDBServer.socket.socket.connect", side_effect=Exception
    )
    def test_default_servers_unavailable(self, _):
        """Check that exception is raised is none of the default server is available."""
        with self.assertRaises(PDBServer.PDBServersConnectionError):
            PDBServer.get_fastest_server(PDBServer.PDBServerProtocol.FTP)


class TestHandleLegacyServer(unittest.TestCase):
    """Test PDBServer.handle_legacy_server."""

    def test_previous_default_value(self):
        """Check that the previous default value of the server params is handled."""
        server, protocol = PDBServer.handle_legacy_server("ftp://ftp.wwpdb.org")
        assert server.domain == "wwpdb.org"
        assert protocol == PDBServer.PDBServerProtocol.FTP

    def test_previous_another_server(self):
        """Check that the previous default value of the server params is handled."""
        server, protocol = PDBServer.handle_legacy_server("https://ebi.ac.uk/")
        assert server.domain == "ebi.ac.uk"
        assert protocol == PDBServer.PDBServerProtocol.HTTPS

    def test_unhandled_protocol(self):
        """Check that exception is raised if the protocol is not handled."""
        with self.assertRaises(PDBServer.UnsupportedProtocolError):
            PDBServer.handle_legacy_server("ftps://ftp.wwpdb.org")

    def test_unhandled_server(self):
        """Check that exception is raised if the server is not handled."""
        with self.assertRaises(PDBServer.UnsupportedServerError):
            PDBServer.handle_legacy_server("ftp://ftp.not-a-supported-server.org")


class RetrievePDBFileTestMixin:
    """Mixin providing methods to test PDB file retrieving using FTP and HTTPS."""

    SERVER: PDBServer
    PDB_CODE = "127d"
    PDB_BUNDLE_CODE = "3k1q"

    def _retrieve_file(self, file_format, protocol, code=None, index=None):
        if (
            self.SERVER not in file_format.servers
            or protocol not in file_format.protocols
        ):
            return
        with tempfile.TemporaryDirectory(dir=os.getcwd()) as temporary_directory_path:
            filepath = pathlib.Path(temporary_directory_path, str(uuid.uuid4()))
            self.SERVER.retrieve_file(
                path=filepath,
                code=code or self.PDB_CODE,
                file_format=file_format,
                protocol=protocol,
                index=index,
            )
            self.assertTrue(
                filepath.exists(),
                msg=f"Error downloading file (format: {file_format}, protocol: {protocol}, server: {self.SERVER}).",
            )

    def test_HTTPS_PDB(self):
        self._retrieve_file(
            protocol=PDBServer.PDBServerProtocol.HTTPS, file_format=FILE_FORMATS["PDB"]
        )

    def test_HTTPS_PDB_ASSEMBLY(self):
        self._retrieve_file(
            protocol=PDBServer.PDBServerProtocol.HTTPS,
            file_format=FILE_FORMATS["PDB_ASSEMBLY"],
            index=1,
        )

    def test_HTTPS_MMCIF(self):
        self._retrieve_file(
            protocol=PDBServer.PDBServerProtocol.HTTPS,
            file_format=FILE_FORMATS["MMCIF"],
        )

    def test_HTTPS_MMCIF_ASSEMBLY(self):
        self._retrieve_file(
            protocol=PDBServer.PDBServerProtocol.HTTPS,
            file_format=FILE_FORMATS["MMCIF_ASSEMBLY"],
            index=1,
        )

    def test_HTTPS_PDBML(self):
        self._retrieve_file(
            protocol=PDBServer.PDBServerProtocol.HTTPS,
            file_format=FILE_FORMATS["PDBML"],
        )

    def test_HTTPS_PDB_BUNDLE(self):
        self._retrieve_file(
            protocol=PDBServer.PDBServerProtocol.HTTPS,
            file_format=FILE_FORMATS["PDB_BUNDLE"],
            code=self.PDB_BUNDLE_CODE,
        )

    def test_HTTPS_MMTF(self):
        self._retrieve_file(
            protocol=PDBServer.PDBServerProtocol.HTTPS, file_format=FILE_FORMATS["MMTF"]
        )

    def test_FTP_PDB(self):
        self._retrieve_file(
            protocol=PDBServer.PDBServerProtocol.FTP, file_format=FILE_FORMATS["PDB"]
        )

    def test_FTP_PDB_ASSEMBLY(self):
        self._retrieve_file(
            protocol=PDBServer.PDBServerProtocol.FTP,
            file_format=FILE_FORMATS["PDB_ASSEMBLY"],
            index=1,
        )

    def test_FTP_MMCIF(self):
        self._retrieve_file(
            protocol=PDBServer.PDBServerProtocol.FTP, file_format=FILE_FORMATS["MMCIF"]
        )

    def test_FTP_MMCIF_ASSEMBLY(self):
        self._retrieve_file(
            protocol=PDBServer.PDBServerProtocol.FTP,
            file_format=FILE_FORMATS["MMCIF_ASSEMBLY"],
            index=1,
        )

    def test_FTP_PDBML(self):
        self._retrieve_file(
            protocol=PDBServer.PDBServerProtocol.FTP, file_format=FILE_FORMATS["PDBML"]
        )

    def test_FTP_PDB_BUNDLE(self):
        self._retrieve_file(
            protocol=PDBServer.PDBServerProtocol.FTP,
            file_format=FILE_FORMATS["PDB_BUNDLE"],
            code=self.PDB_BUNDLE_CODE,
        )

    def test_FTP_MMTF(self):
        self._retrieve_file(
            protocol=PDBServer.PDBServerProtocol.FTP, file_format=FILE_FORMATS["MMTF"]
        )


class TestDownloadPDBFromWWServer(RetrievePDBFileTestMixin, unittest.TestCase):
    """Test PDB file retrieving using FTP and HTTPS on the worldwide server."""

    SERVER = PDBServer.PDB_SERVERS["WW"]


class TestDownloadPDBFromUSServer(RetrievePDBFileTestMixin, unittest.TestCase):
    """Test PDB file retrieving using FTP and HTTPS on the US server."""

    SERVER = PDBServer.PDB_SERVERS["US"]


class TestDownloadPDBFromUKServer(RetrievePDBFileTestMixin, unittest.TestCase):
    """Test PDB file retrieving using FTP and HTTPS on the UK server."""

    SERVER = PDBServer.PDB_SERVERS["UK"]


class TestDownloadPDBFromJPServer(RetrievePDBFileTestMixin, unittest.TestCase):
    """Test PDB file retrieving using FTP and HTTPS on the JP server."""

    SERVER = PDBServer.PDB_SERVERS["JP"]


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
