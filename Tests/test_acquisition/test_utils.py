"""Test acquisition.utils."""
import unittest
import unittest.mock

from Bio.PDB.acquisition.utils import build_socket


class TestAcquisitionUtilsBuildSocket(unittest.TestCase):
    """Test acquisition.utils.build_socket."""

    def test_build_socket_auto_close(self):
        """Check that the socket is closed exiting the with."""
        with unittest.mock.patch(
            "Bio.PDB.acquisition.utils.socket.socket.close"
        ) as mock_close_socket:
            with build_socket() as server_socket:
                pass
            mock_close_socket.assert_called_once_with()

    def test_build_socket_auto_close_on_exception(self):
        """Check that the socket is closed exiting the with even if an exception occurs."""
        with unittest.mock.patch(
            "Bio.PDB.acquisition.utils.socket.socket.close"
        ) as mock_close_socket:
            try:
                with build_socket() as server_socket:
                    raise Exception
            except Exception:
                pass
            mock_close_socket.assert_called_once_with()


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
