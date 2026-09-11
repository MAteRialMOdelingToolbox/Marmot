import argparse
import os

import numpy as np


def _extract_strain(history, strain_field=None):
    if strain_field is not None:
        return np.array([getattr(h, strain_field) for h in history])
    if hasattr(history[0], "F"):
        return np.array([h.F for h in history])
    if hasattr(history[0], "strain"):
        return np.array([h.strain for h in history])
    raise AttributeError("HistoryEntry has neither 'F' nor 'strain' attribute.")


def writeReferenceSolution(history, reference_file, strain_field=None):
    """Write solver history to a reference npz file."""
    if len(history) == 0:
        raise RuntimeError("Solver history is empty!")

    stress = np.array([h.stress for h in history])
    strain = _extract_strain(history, strain_field)

    np.savez(reference_file, stress=stress, strain=strain)
    print(f"Reference solution created: {reference_file}")


def compareHistoryToReferenceSolution(history, reference_file, strain_field=None):
    """Compare a solver history against a reference npz file."""
    if len(history) == 0:
        raise RuntimeError("Solver history is empty!")

    if not os.path.exists(reference_file):
        raise FileNotFoundError(
            f"Reference file not found: {reference_file}. Run with --writeReferenceSolution to create it."
        )

    stress = np.array([h.stress for h in history])
    strain = _extract_strain(history, strain_field)

    ref = np.load(reference_file)
    np.testing.assert_allclose(stress, ref["stress"], atol=1e-6, rtol=1e-5)
    np.testing.assert_allclose(strain, ref["strain"], atol=1e-6, rtol=1e-5)
    print(f"Reference solution passed: {reference_file}")


def run_test(history, script_file, strain_field=None):
    """
    Unified entry point for testing python scripts.
    Parses command line arguments to either compare or write the reference solution.
    """
    parser = argparse.ArgumentParser(description="Marmot example test script")
    parser.add_argument("--writeReferenceSolution", action="store_true", help="Write the reference solution")

    args, _ = parser.parse_known_args()

    reference_file = os.path.splitext(script_file)[0] + "_reference.npz"

    if args.writeReferenceSolution:
        writeReferenceSolution(history, reference_file, strain_field=strain_field)
    else:
        compareHistoryToReferenceSolution(history, reference_file, strain_field=strain_field)
