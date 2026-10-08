#!/usr/bin/env bash
# Settings and input list are in combine_edr_instruction.txt beside this script.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

fail() { printf 'Error: %s\n' "$*" >&2; exit 1; }
DRY_RUN=false
case "${1:-}" in
    '') ;;
    --dry-run) DRY_RUN=true ;;
    *) fail "Usage: bash combine_edr.sh [--dry-run]" ;;
esac
[[ $# -le 1 ]] || fail "Too many arguments."
[[ -f combine_edr_instruction.txt ]] || fail "Missing combine_edr_instruction.txt."

GMX=gmx
OUTPUT=
FIRST=
PARTS_PREFIX=
inputs=()
while IFS= read -r line || [[ -n "$line" ]]; do
    line="${line%$'\r'}"
    [[ "$line" =~ ^[[:space:]]*$ || "$line" =~ ^[[:space:]]*# ]] && continue
    [[ "$line" == *=* ]] || fail "Expected KEY=value: $line"
    key="${line%%=*}"
    value="${line#*=}"
    [[ -n "$value" ]] || fail "Empty setting: $key"
    case "$key" in
        GMX) GMX="$value" ;;
        OUTPUT) OUTPUT="$value" ;;
        FIRST) FIRST="$value" ;;
        PARTS_PREFIX) PARTS_PREFIX="$value" ;;
        INPUT) inputs+=("$value") ;;
        *) fail "Unknown setting: $key" ;;
    esac
done < combine_edr_instruction.txt

if [[ -n "$FIRST" || -n "$PARTS_PREFIX" ]]; then
    [[ -n "$FIRST" && -n "$PARTS_PREFIX" ]] || fail "Set both FIRST and PARTS_PREFIX."
    [[ ${#inputs[@]} -eq 0 ]] || fail "Use FIRST/PARTS_PREFIX or manual INPUT lines, not both."
    # Expand only this production prefix and four-digit part numbers.
    # Missing trailing parts are fine; unmatched patterns create no filenames.
    export LC_ALL=C
    shopt -s nullglob
    parts=("${PARTS_PREFIX}".part[0-9][0-9][0-9][0-9].edr)
    [[ ${#parts[@]} -gt 0 ]] || fail "No numbered EDR parts found for $PARTS_PREFIX."
    inputs=("$FIRST" "${parts[@]}")
fi
[[ ${#inputs[@]} -ge 2 ]] || fail "At least two input files are needed."
[[ "$OUTPUT" == *.edr ]] || fail "OUTPUT must end in .edr."
for file in "${inputs[@]}"; do
    [[ "$file" == *.edr && -f "$file" && -r "$file" && -s "$file" ]] || fail "Missing, empty or unreadable EDR: $file"
done

command=("$GMX" eneconv -f "${inputs[@]}" -o "$OUTPUT")
printf 'Found %s input EDR files:\n' "${#inputs[@]}"
printf '  %s\n' "${inputs[@]}"
if "$DRY_RUN"; then
    printf '%q ' "${command[@]}"; printf '\n'
    printf 'Dry run only: no files written.\n'
    if [[ -e "$OUTPUT" || -L "$OUTPUT" ]]; then
        printf 'Output already exists; a real run will stop: %s\n' "$OUTPUT"
    fi
    exit 0
fi
[[ ! -e "$OUTPUT" && ! -L "$OUTPUT" ]] || fail "Output already exists: $OUTPUT. Choose a new OUTPUT filename."
command -v "$GMX" >/dev/null || fail "Cannot find $GMX. Activate your GROMACS environment."
"${command[@]}"
[[ -s "$OUTPUT" ]] || fail "No non-empty output was produced. Inspect any partial output before retrying."
printf 'Combined EDR saved to %s\n' "$OUTPUT"
