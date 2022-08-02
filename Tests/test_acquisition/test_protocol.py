"""Test acquisition.protocol."""
import unittest
import unittest.mock

from Bio.PDB.acquisition import exceptions
from Bio.PDB.acquisition.protocol import PDBServerProtocol, get_default_protocol


class TestPDBServerProtocol(unittest.TestCase):
    """Test acquisition.protocol.PDBServerProtocol."""

    def test_url_prefix(self):
        """Test PDBServerProtocol.url_prefix."""
        self.assertEqual(PDBServerProtocol.FTP.url_prefix, "ftp://")
        self.assertEqual(PDBServerProtocol.HTTPS.url_prefix, "https://")


class TestPDBServerProtocolFunctions(unittest.TestCase):
    """Test acquisition.protocol module functions."""

    def test_get_default_protocol(self):
        """Test protocol.get_default_protocol."""
        config = {"default": {"protocol": "HTTPS"}}
        with unittest.mock.patch(
            "Bio.PDB.acquisition.protocol.get_config", return_value=config
        ):
            self.assertEqual(get_default_protocol(), PDBServerProtocol.HTTPS)

    def test_get_default_protocol_error(self):
        """Test protocol.get_default_protocol."""
        missing_default_protocol_config = {"default": {}}
        with (
            unittest.mock.patch(
                "Bio.PDB.acquisition.protocol.get_config",
                return_value=missing_default_protocol_config,
            ),
            self.assertRaises(exceptions.ConfigurationError),
        ):
            get_default_protocol()


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
