import argparse
import subprocess
import sys
from pathlib import Path

parser = argparse.ArgumentParser(description="Build the Marmot documentation (Doxygen + Sphinx).")
parser.add_argument(
    "--werror",
    action="store_true",
    help="fail if Doxygen or Sphinx report any warning (as in CI)",
)
args = parser.parse_args()

# change to doc directory
doc_dir = Path(__file__).resolve().parent.parent / "doc"

# run doxygen; with --werror, the Doxyfile is passed on stdin with WARN_AS_ERROR overridden
if args.werror:
    doxyfile = (doc_dir / "Doxyfile").read_text() + "\nWARN_AS_ERROR = FAIL_ON_WARNINGS\n"
    subprocess.run(["doxygen", "-"], cwd=doc_dir, input=doxyfile, text=True, check=True)
else:
    subprocess.run(["doxygen", "Doxyfile"], cwd=doc_dir, check=True)

# run sphinx; with --werror, warnings are errors, but all of them are reported before failing
subprocess.run(
    [
        sys.executable,
        "-m",
        "sphinx",
        "-b",
        "html",
        *(["-W", "--keep-going"] if args.werror else []),
        "-Dbreathe_projects.Marmot=doc_out/xml",
        ".",
        "doc_out/sphinx",
    ],
    cwd=doc_dir,
    check=True,
)
