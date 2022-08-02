"""Test acquisition.legacy."""
import unittest
import unittest.mock

from Bio.PDB.acquisition import exceptions
from Bio.PDB.acquisition.legacy import handle_legacy_server
from Bio.PDB.acquisition.protocol import PDBServerProtocol


class TestAcquisitionLegacy(unittest.TestCase):
    """Test acquisition.legacy."""

    def test_previous_default_value(self):
        """Check that the previous default value of the server params is handled."""
        server, protocol = handle_legacy_server("ftp://ftp.wwpdb.org")
        assert server.domain == "wwpdb.org"
        assert protocol == PDBServerProtocol.FTP

    def test_previous_another_server(self):
        """Check that the previous default value of the server params is handled."""
        server, protocol = handle_legacy_server("https://ebi.ac.uk/")
        assert server.domain == "ebi.ac.uk"
        assert protocol == PDBServerProtocol.HTTPS

    def test_unhandled_protocol(self):
        """Check that exception is raised if the protocol is not handled."""
        with self.assertRaises(exceptions.UnsupportedProtocolError):
            handle_legacy_server("ftps://ftp.wwpdb.org")

    def test_unhandled_server(self):
        """Check that exception is raised if the server is not handled."""
        with self.assertRaises(exceptions.UnsupportedServerError):
            handle_legacy_server("ftp://ftp.not-a-supported-server.org")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
