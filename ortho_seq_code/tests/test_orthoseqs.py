import numpy as np
import os

from click.testing import CliRunner

from ortho_seq_code.orthogonal_polynomial import orthogonal_polynomial
from ortho_seq_code.cli import cli
from ortho_seq_code.tests import orthoseqs_tst_utils as utils


def test_cli(protein_seqs_no_padding, protein_pheno_no_padding):
    molecule = "protein"
    sites = 6
    dm = 20
    pop_size = 6
    poly_order = "first"
    out_dir = "/tmp"

    runner = CliRunner()

    result = runner.invoke(
        cli,
        [
            protein_seqs_no_padding,
            "--pheno_file",
            protein_pheno_no_padding,
            "--molecule",
            molecule,
            "--sites",
            sites,
            "--dm",
            dm,
            "--pop_size",
            pop_size,
            "--poly_order",
            poly_order,
            "--out_dir",
            out_dir,
        ],
    )

    assert result.exit_code == 0


def assert_equality(expected_path, actual_path):
    assert os.path.exists(expected_path)
    assert os.path.exists(actual_path)
    obtained_arrays = np.load(actual_path)
    for key, obtained_array in obtained_arrays.items():
        with np.load(expected_path) as expected_arrays:
            expected_array = expected_arrays[key]
            np.testing.assert_array_equal(expected_array, obtained_array)


def test_nucleotide_first_order(
        nucleotide_first_order_data_dir, nucleotide_params_first_order):

    with utils.TempDirectory() as location:
        nucleotide_params_first_order = nucleotide_params_first_order._replace(
            out_dir=location
        )
        orthogonal_polynomial(*nucleotide_params_first_order)
        basename = os.path.basename(nucleotide_params_first_order.seqs_filename)
        expected_path = os.path.join(nucleotide_first_order_data_dir, basename + ".npz")
        obtained_path = os.path.join(location, basename + ".npz")
        assert_equality(expected_path, obtained_path)


def test_nucleotide_second_order(
        nucleotide_second_order_data_dir, nucleotide_params_second_order):
    with utils.TempDirectory() as location:
        nucleotide_params_second_order = nucleotide_params_second_order._replace(
            out_dir=location
        )
        orthogonal_polynomial(*nucleotide_params_second_order)

        basename = os.path.basename(nucleotide_params_second_order.seqs_filename)
        expected_path = os.path.join(
            nucleotide_second_order_data_dir, basename + ".npz"
        )
        obtained_path = os.path.join(location, basename + ".npz")
        assert_equality(expected_path, obtained_path)


def test_protein_first_order(protein_data_dir, protein_params_first_order):

    with utils.TempDirectory() as location:
        protein_params_first_order = protein_params_first_order._replace(
            out_dir=location
        )
        orthogonal_polynomial(*protein_params_first_order)

        basename = os.path.basename(protein_params_first_order.seqs_filename)
        expected_path = os.path.join(protein_data_dir, basename + ".npz")
        obtained_path = os.path.join(location, basename + ".npz")
        assert_equality(expected_path, obtained_path)
