# Copyright (C) 2022, Julian Maurin (julian.maurin.dev@proton.me)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Testing PDB server utils."""

import unittest
import unittest.mock

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
    @unittest.mock.patch("Bio.PDB.PDBServer.PDB_SERVERS", PDBServer.PDB_SERVERS[:2])
    @unittest.mock.patch("Bio.PDB.PDBServer.socket.socket.connect", return_value=None)
    def test_timing(self, _, __):
        """Check that the returned server is the fastest.

        Mocking the first server duration to 5s (5-0) and the second to 1s (6-5).
        """
        self.assertEqual(
            PDBServer.get_fastest_server(PDBServer.PDBServerProtocol.FTP),
            PDBServer.PDB_SERVERS[1],
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


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
