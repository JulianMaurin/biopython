"""Test acquisition.server.PDBServer."""
import dataclasses
import unittest
import unittest.mock

from Bio.PDB import acquisition
from Bio.PDB.acquisition import exceptions
from Bio.PDB.acquisition.protocol import PDBServerProtocol
from Bio.PDB.acquisition.server import _get_server_connection_timing, get_pdb_servers

WWServer = acquisition.PDBServer(
    code="WW", label="Worldwide", domain="wwpdb.org", pdb_directory="/pub/pdb"
)
USServer = acquisition.PDBServer(
    code="US", label="United States", domain="rcsb.org", pdb_directory="/pub/pdb"
)


class TestPDBServer(unittest.TestCase):
    """Test acquisition.server.PDBServer."""

    def test_immutable(self):
        """Test PDBServer dataclass immutability."""
        with self.assertRaises(dataclasses.FrozenInstanceError):
            WWServer.domain = "test"

    def test_build_pdb_directory_url(self):
        """Test PDBServer.build_pdb_directory_url."""
        self.assertEqual(
            WWServer.build_pdb_directory_url(protocol=PDBServerProtocol.FTP),
            "ftp://ftp.wwpdb.org/pub/pdb",
        )

    def test_subdomain(self):
        """Test PDBServer.subdomain."""
        self.assertEqual(
            acquisition.PDBServer.subdomain(protocol=PDBServerProtocol.FTP), "ftp."
        )
        self.assertEqual(
            acquisition.PDBServer.subdomain(protocol=PDBServerProtocol.HTTPS), ""
        )


class TestPDBServerFunctions(unittest.TestCase):
    """Test acquisition.server module functions."""

    def test_get_pdb_servers(self):
        """Test server.get_pdb_servers."""
        config = {"servers": {"A": {"label": "B", "domain": "C", "pdb_directory": "D"}}}
        with unittest.mock.patch(
            "Bio.PDB.acquisition.server.get_config", return_value=config
        ):
            servers = get_pdb_servers()
        self.assertEqual(len(servers), 1)
        self.assertIsInstance(servers[0], acquisition.PDBServer)
        self.assertEqual(servers[0].code, "A")
        self.assertEqual(servers[0].label, "B")
        self.assertEqual(servers[0].domain, "C")
        self.assertEqual(servers[0].pdb_directory, "D")

    def test_get_pdb_servers_config_error(self):
        missing_domain_config = {"servers": {"A": {"label": "B", "pdb_directory": "D"}}}
        with (
            unittest.mock.patch(
                "Bio.PDB.acquisition.server.get_config",
                return_value=missing_domain_config,
            ),
            self.assertRaises(exceptions.ConfigurationError),
        ):
            get_pdb_servers()

    @unittest.mock.patch(
        "Bio.PDB.acquisition.utils.socket.socket.connect", unittest.mock.MagicMock()
    )
    def test_get_server_connection_timing(self):
        """Test server._get_server_connection_timing."""
        _get_server_connection_timing.cache_clear()
        with (
            unittest.mock.patch(
                "Bio.PDB.acquisition.server.time.time", side_effect=[2, 5, 5, 6]
            ),
        ):
            # 5s - 2s, the timing is equals to 3s
            self.assertEqual(
                _get_server_connection_timing(WWServer, PDBServerProtocol.FTP), 3
            )
            # the result for these parameters is stored in the cache, time is not called
            self.assertEqual(
                _get_server_connection_timing(WWServer, PDBServerProtocol.FTP), 3
            )
            # 6s - 5s, the timing is equals to 1s
            self.assertEqual(
                _get_server_connection_timing(WWServer, PDBServerProtocol.HTTPS), 1
            )

    def test_get_fastest_server(self):
        with (
            unittest.mock.patch(
                "Bio.PDB.acquisition.server.get_pdb_servers",
                return_value=[WWServer, USServer],
            ),
            unittest.mock.patch(
                "Bio.PDB.acquisition.server._get_server_connection_timing",
                side_effect=[5, 1],
            ),
        ):
            self.assertEqual(
                acquisition.get_fastest_server(PDBServerProtocol.FTP), USServer
            )

    @unittest.mock.patch(
        "Bio.PDB.acquisition.utils.socket.socket.connect", side_effect=Exception
    )
    def test_default_servers_unavailable(self, _):
        """Check that exception is raised is none of the default server is available."""
        with self.assertRaises(exceptions.PDBServersConnectionError):
            acquisition.get_fastest_server(PDBServerProtocol.FTP)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
