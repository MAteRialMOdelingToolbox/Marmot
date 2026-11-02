"""Patch nanobind so that zero-argument functions do not declare a zero-sized array.

nanobind declares ``arg_data args[Size];`` in ``nb_attr.h``. For a function without
arguments ``Size`` is 0, and a zero-sized array is not valid ISO C++, which makes
compilers reject the header in pedantic mode. The patch clamps the extent to 1.

Called as a ``PATCH_COMMAND`` from FetchContent, so it must be idempotent: a source
tree that has already been patched is left untouched. If neither the original nor the
patched pattern is present, the vendored nanobind version has changed and the script
fails, rather than silently producing an unpatched build.
"""

import sys

ORIGINAL = "    // GCC and Clang do.\n    arg_data args[Size];"
PATCHED = "    // GCC and Clang do.\n    arg_data args[Size == 0 ? 1 : Size];"

path = sys.argv[1]

with open(path, "r") as f:
    content = f.read()

if ORIGINAL not in content:
    if PATCHED in content:
        print(f"{path} is already patched, nothing to do.")
        sys.exit(0)
    print(
        f"error: expected pattern '{ORIGINAL}' not found in {path}. "
        "The nanobind version has probably changed; patch_nanobind.py needs to be updated.",
        file=sys.stderr,
    )
    sys.exit(1)

with open(path, "w") as f:
    f.write(content.replace(ORIGINAL, PATCHED))

print(f"patched {path}")
