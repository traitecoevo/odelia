#!/usr/bin/env bash
#
# Re-apply odelia's local patches to the vendored XAD library.
#
# The XAD sources under inst/include/XAD/ and src/Tape.cpp are vendored from
# https://github.com/auto-differentiation/xad. We carry a small number of local
# patches on top of that snapshot (see tools/patches/ and tools/README.md).
#
# Run this after re-vendoring XAD so the local fixes are restored:
#
#   tools/apply-xad-patches.sh          # apply all patches
#   tools/apply-xad-patches.sh --check  # dry-run; verify they still apply
#
# The script is safe to re-run: patches that are already applied are skipped.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
patch_dir="${repo_root}/tools/patches"

cd "${repo_root}"

check_only=0
if [[ "${1:-}" == "--check" ]]; then
  check_only=1
fi

shopt -s nullglob
patches=("${patch_dir}"/*.patch)
if [[ ${#patches[@]} -eq 0 ]]; then
  echo "No patches found in ${patch_dir}"
  exit 0
fi

for p in "${patches[@]}"; do
  name="$(basename "${p}")"
  if git apply --reverse --check "${p}" >/dev/null 2>&1; then
    echo "skip   ${name} (already applied)"
    continue
  fi
  if ! git apply --check "${p}" >/dev/null 2>&1; then
    echo "FAIL   ${name} (does not apply cleanly - upstream likely changed; re-create the patch)" >&2
    exit 1
  fi
  if [[ ${check_only} -eq 1 ]]; then
    echo "ok     ${name} (would apply cleanly)"
  else
    git apply "${p}"
    echo "apply  ${name}"
  fi
done

echo "Done."
