import unittest
from unittest.mock import mock_open, patch, MagicMock
from pathlib import Path
import requests

import numpy as np

from kinex.functions import (
    check_sequence,
    get_sequence_format,
    download_file_to_resource,
    get_distances,
)


class TestIO(unittest.TestCase):

    def test_get_sequence_format(self):
        self.assertEqual(get_sequence_format("PSVEPPLs*QETFSDL"), "*")
        self.assertEqual(get_sequence_format("PSVEXPLs*QXTF___"), "*")
        self.assertEqual(get_sequence_format("PSVEPPLsQETFSDL"), "central")
        self.assertEqual(get_sequence_format("PSVEPPLs(ph)QETFSDL"), "(ph)")
        self.assertEqual(get_sequence_format("LQVKIPSKEEEsAD"), "unsupported")

    def test_check_sequence(self):
        self.assertEqual(check_sequence("PSVEPPLs*QETFSDL", sequence_format="*"), True)
        self.assertEqual(check_sequence("PSVEXPLs*QXTF___", sequence_format="*"), True)
        self.assertEqual(
            check_sequence("PSVEPPLsQETFSDL", sequence_format="central"), True
        )

    def test_valid_input(self):
        experiment1 = {
            "dominant_enrichment_value_log2": [1, 2, 3],
            "dominant_p_value_log10_abs": [4, 5, 6],
        }
        experiment2 = {
            "dominant_enrichment_value_log2": [1, 1, 1],
            "dominant_p_value_log10_abs": [4, 4, 4],
        }
        expected = np.sqrt(np.array([0, 1, 4]) + np.array([0, 1, 4]))
        result = get_distances(experiment1, experiment2)
        np.testing.assert_array_almost_equal(result, expected)

    def test_different_shapes(self):
        experiment1 = {
            "dominant_enrichment_value_log2": [1, 2],
            "dominant_p_value_log10_abs": [4, 5],
        }
        experiment2 = {
            "dominant_enrichment_value_log2": [1, 1, 1],
            "dominant_p_value_log10_abs": [4, 4, 4],
        }
        with self.assertRaises(ValueError):
            get_distances(experiment1, experiment2)

    def test_missing_key(self):
        experiment1 = {"dominant_enrichment_value_log2": [1, 2, 3]}
        experiment2 = {
            "dominant_enrichment_value_log2": [1, 1, 1],
            "dominant_p_value_log10_abs": [4, 4, 4],
        }
        with self.assertRaises(ValueError):
            get_distances(experiment1, experiment2)

    def test_empty_input(self):
        experiment1 = {
            "dominant_enrichment_value_log2": [],
            "dominant_p_value_log10_abs": [],
        }
        experiment2 = {
            "dominant_enrichment_value_log2": [],
            "dominant_p_value_log10_abs": [],
        }
        result = get_distances(experiment1, experiment2)
        np.testing.assert_array_equal(result, np.array([]))

    @patch("requests.get")
    @patch("kinex.functions.resources.path")
    @patch("builtins.open", new_callable=mock_open)
    def test_download_file_to_resource_success(
        self, mock_open, mock_resources_path, mock_requests_get
    ):
        # Mock the path returned by resources.path as a Path object
        mock_file_path = MagicMock()
        mock_file_path.__enter__.return_value = Path("/mocked/path/to/resource.csv.gz")
        mock_resources_path.return_value = mock_file_path

        # Mock the requests.get response
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.content = b"mocked file content"
        mock_requests_get.return_value = mock_response

        # Call the function to test
        download_file_to_resource("http://example.com/file.csv.gz", "resource.csv.gz")

        # Assertions
        mock_requests_get.assert_called_once_with(
            "http://example.com/file.csv.gz", stream=True, timeout=10
        )
        mock_resources_path.assert_called_once_with(
            "kinex.resources", "resource.csv.gz"
        )

        # Verify the open() method was called to write the content
        mock_open.assert_called_once_with(Path("/mocked/path/to/resource.csv.gz"), "wb")
        mock_open().write.assert_called_once_with(b"mocked file content")

    # Download Failure Due to Network Issues
    @patch("requests.get")
    @patch("kinex.functions.resources.path")
    @patch("builtins.open", new_callable=mock_open)
    def test_download_file_to_resource_network_failure(
        self, mock_open, mock_resources_path, mock_requests_get
    ):
        # Mock the path returned by resources.path
        mock_file_path = MagicMock()
        mock_file_path.__enter__.return_value = Path("/mocked/path/to/resource.csv.gz")
        mock_resources_path.return_value = mock_file_path

        # Simulate a network error during requests.get
        mock_requests_get.side_effect = requests.exceptions.ConnectionError(
            "Network error"
        )

        # Assert that the function raises an exception
        with self.assertRaises(requests.exceptions.ConnectionError):
            download_file_to_resource(
                "http://example.com/file.csv.gz", "resource.csv.gz"
            )

        # Ensure no file operations occurred
        mock_open.assert_not_called()
        mock_resources_path.assert_called_once_with(
            "kinex.resources", "resource.csv.gz"
        )

    # Test the case where Invalid URL is passed
    @patch("requests.get")
    @patch("kinex.functions.resources.path")
    @patch("builtins.open", new_callable=mock_open)
    def test_download_file_to_resource_invalid_url(
        self, mock_open, mock_resources_path, mock_requests_get
    ):
        # Mock the path returned by resources.path
        mock_file_path = MagicMock()
        mock_file_path.__enter__.return_value = Path("/mocked/path/to/resource.csv.gz")
        mock_resources_path.return_value = mock_file_path

        # Simulate a requests.exceptions.MissingSchema exception
        mock_requests_get.side_effect = requests.exceptions.MissingSchema("Invalid URL")

        # Assert that a ValueError is raised
        with self.assertRaises(ValueError):
            download_file_to_resource("invalid_url", "resource.csv.gz")

        # Ensure no file operations occurred
        mock_open.assert_not_called()
        mock_resources_path.assert_called_once_with(
            "kinex.resources", "resource.csv.gz"
        )


if __name__ == "__main__":
    unittest.main()
