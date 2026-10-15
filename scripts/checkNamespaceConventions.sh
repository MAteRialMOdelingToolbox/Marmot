#!/usr/bin/env bash
# Check the Marmot namespace convention (see AGENTS.md, "Namespace Convention").
# Usage: scripts/checkNamespaceConventions.sh [--warn]
#   --warn  report violations but exit 0
set -u
cd "$(dirname "$0")/.."

MODE_EXIT=1
[ "${1:-}" = "--warn" ] && MODE_EXIT=0

FAIL=0
report() { # title, matches
  if [ -n "$2" ]; then
    echo "== $1"
    echo "$2"
    FAIL=1
  fi
}

HDR=$(find modules -name "*.h" -not -path "*/build*")
ALL=$(find modules tests -name "*.h" -o -name "*.cpp" | grep -v build)

# 1. The MarmotLibrary namespace must not reappear.
report "MarmotLibrary namespace is abolished (use Marmot::Registry)" \
  "$(grep -ln "MarmotLibrary" $ALL 2>/dev/null)"

# 2. No classes/structs in the global namespace of installed headers
#    (heuristic: class/struct at column 0; forward declarations excluded).
report "class/struct in global namespace in a header" \
  "$(grep -lnE "^(class|struct) +(MARMOT_API +)?[A-Za-z_]+[^;]*$" $HDR 2>/dev/null)"

# 3. Old derivative-order spelling.
report "use FirstOrder/SecondOrder/ThirdOrder, not *OrderDerived namespaces" \
  "$(grep -lnE "namespace (First|Second|Third)OrderDerived|OrderDerived::" $ALL 2>/dev/null)"

# 4. Lowercase namespace names (PascalCase required).
report "lowercase namespace name" \
  "$(grep -nE "^\s*namespace +[a-z][a-zA-Z0-9_]* *\{" $ALL 2>/dev/null | grep -v "namespace *{")"

# 5. using-directives / using-declarations at global scope in headers
#    (column 0 => global scope; function/class scope is indented).
report "using at global scope in a header" \
  "$(grep -nE "^using " $HDR 2>/dev/null)"

# 6. Old-style nested namespace declarations (C++17 required): a bare
#    'namespace Marmot {' opener whose children are opened separately.
#    (Sibling sub-namespace groups inside one qualified opener are fine.)
report "old-style 'namespace Marmot {' opener (use C++17 'namespace Marmot::X {')" \
  "$(grep -lnE "^namespace Marmot \{" $ALL 2>/dev/null)"

# 7. Namespace depth > 6 (counting Marmot).
report "namespace deeper than 6 levels" \
  "$(grep -nE "namespace +(\w+::){6,}\w+ *\{" $ALL 2>/dev/null)"

if [ $FAIL -ne 0 ]; then
  echo "Namespace convention violations found."
  exit $MODE_EXIT
fi
echo "Namespace conventions OK."
