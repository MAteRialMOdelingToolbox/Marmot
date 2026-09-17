import subprocess
import sys
from pathlib import Path

# change to doc directory
doc_dir = Path(__file__).resolve().parent.parent / "doc"

# run doxygen
subprocess.run(["doxygen", "Doxyfile"], cwd=doc_dir, check=True)

# run sphinx
subprocess.run(
    [
        sys.executable,
        "-m",
        "sphinx",
        "-b",
        "html",
        "-Dbreathe_projects.Marmot=doc_out/xml",
        ".",
        "doc_out/sphinx",
    ],
    cwd=doc_dir,
    check=True,
)
