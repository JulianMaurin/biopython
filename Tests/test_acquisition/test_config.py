"""Test acquisition.config."""
import unittest
import unittest.mock
from Bio.PDB.acquisition.config import get_config


class TestAcquisitionConfig(unittest.TestCase):
    """Test acquisition.config."""

    def test_get_config(self):
        """Test config.get_config."""
        config = get_config()
        self.assertIn("default", config)
        self.assertIn("protocol", config["default"])
        self.assertIn("servers", config)
        self.assertEqual(len(config["servers"]), 4)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
