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
ALL=$(find modules tests examples -name "*.h" -o -name "*.cpp" | grep -v build)

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

# 4. Lowercase namespace names (PascalCase required), checking every component
#    of a C++17 nested opener. Anonymous namespaces never match (an identifier
#    is required after 'namespace').
report "lowercase namespace name" \
  "$(grep -nE "^\s*namespace +[A-Za-z_][A-Za-z0-9_]*(::[A-Za-z_][A-Za-z0-9_]*)* *\{" $ALL 2>/dev/null \
    | grep -E "namespace +[a-z]|::[a-z]" \
    | grep -vE "namespace +(nanobind|std|autodiff|nb)\b")"

# 5. using-directives / using-declarations at global or namespace scope in
#    headers. Tracks braces: a 'using' line (aliases 'using X = ...' excluded)
#    is flagged unless some enclosing brace is a non-namespace block (class,
#    function, ...). Comment-only lines are skipped so braces in docs do not
#    skew the depth. AGENTS.md rule 9 sanctions the autodiff import in
#    Marmot::Math::Differentiation::Automatic (the MarmotAutomaticDifferentiation*
#    headers).
report "using at global or namespace scope in a header" \
  "$(awk '
      FNR == 1 { depth = 0; nsCount = 0 }
      /^[[:space:]]*(\/\/|\/?\*)/ { next }
      {
        if ($0 ~ /^[[:space:]]*using[[:space:]]/ && $0 !~ /=/ && depth == nsCount)
          print FILENAME ":" FNR ":" $0
        isNs = ($0 ~ /^[[:space:]]*(inline[[:space:]]+)?namespace[^;{]*\{/) ? 1 : 0
        n = length($0)
        for (i = 1; i <= n; i++) {
          c = substr($0, i, 1)
          if (c == "{")      { depth++; ns[depth] = isNs; nsCount += isNs }
          else if (c == "}") { if (depth > 0) { nsCount -= ns[depth]; depth-- } }
        }
      }
    ' $HDR 2>/dev/null | grep -v "MarmotAutomaticDifferentiation")"

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
