#!/usr/bin/env bash
# Read settings from instrcution.txt beside this script; requires only Bash and GROMACS.
set -euo pipefail
script_dir="${BASH_SOURCE[0]%/*}"
[[ "$script_dir" != "${BASH_SOURCE[0]}" ]] || script_dir=.
cd -- "$script_dir"

fail() { printf 'Error: %s\n' "$*" >&2; exit 1; }
DRY_RUN=false
case "${1:-}" in
    '') ;;
    --dry-run) DRY_RUN=true ;;
    *) fail "Usage: bash process_trajectory.sh [--dry-run]" ;;
esac
[[ $# -le 1 ]] || fail "Too many arguments."
[[ -f instrcution.txt ]] || fail "Missing instrcution.txt beside this script."

# Plain KEY=value settings, read as data (never executed as shell commands).
INPUT= TPR= INDEX= CENTER_GROUP= OUTPUT_GROUP= STRIDE= PROCESSED= PREVIEW= GMX=
seen='|'
while IFS= read -r line || [[ -n "$line" ]]; do
    line="${line%$'\r'}"
    [[ "$line" =~ ^[[:space:]]*$ || "$line" =~ ^[[:space:]]*# ]] && continue
    [[ "$line" == *=* ]] || fail "Expected KEY=value in instrcution.txt: $line"
    key="${line%%=*}"
    value="${line#*=}"
    [[ "$seen" != *"|$key|"* ]] || fail "Duplicate setting: $key"
    seen+="$key|"
    [[ -n "$value" ]] || fail "Empty setting: $key"
    case "$key" in
        INPUT|TPR|INDEX|CENTER_GROUP|OUTPUT_GROUP|STRIDE|PROCESSED|PREVIEW|GMX)
            printf -v "$key" '%s' "$value" ;;
        *) fail "Unknown setting: $key" ;;
    esac
done < instrcution.txt

for key in INPUT TPR INDEX CENTER_GROUP OUTPUT_GROUP STRIDE PROCESSED PREVIEW GMX; do
    [[ -n "${!key}" ]] || fail "Missing setting: $key"
done
[[ "$STRIDE" =~ ^[1-9][0-9]*$ ]] || fail "STRIDE must be a positive integer."
command -v "$GMX" >/dev/null || fail "Cannot find $GMX. Activate your GROMACS environment."
for file in "$INPUT" "$TPR" "$INDEX"; do
    [[ -r "$file" && -f "$file" ]] || fail "Missing or unreadable input: $file"
done
[[ "$TPR" == *.tpr ]] || fail "Use a matching production .tpr for -pbc mol."
[[ "$PROCESSED" != "$PREVIEW" ]] || fail "The two output filenames must differ."
for file in "$PROCESSED" "$PREVIEW"; do
    [[ "$file" == *.xtc && "$file" != */* ]] || fail "Outputs must be .xtc filenames inside this folder."
    [[ ! -e "$file" && ! -L "$file" ]] || fail "Output already exists: $file. Rename it or change the output setting."
done

process=("$GMX" trjconv -f "$INPUT" -s "$TPR" -n "$INDEX"
         -pbc mol -ur compact -center -o "$PROCESSED")
preview=("$GMX" trjconv -f "$PROCESSED" -s "$TPR" -n "$INDEX"
         -skip "$STRIDE" -o "$PREVIEW")

printf 'Full-resolution analysis trajectory: %s\nVMD preview: %s\n' "$PROCESSED" "$PREVIEW"
if "$DRY_RUN"; then
    printf '\nStep 1 selections: %s, then %s\n' "$CENTER_GROUP" "$OUTPUT_GROUP"
    printf '%q ' "${process[@]}"; printf '\n'
    printf '\nStep 2 selection: %s\n' "$OUTPUT_GROUP"
    printf '%q ' "${preview[@]}"; printf '\n'
    printf '\nDry run only: no GROMACS commands executed or outputs written.\n'
    exit 0
fi

printf '\nStep 1: centre and compact-wrap all frames...\n'
printf '%s\n%s\n' "$CENTER_GROUP" "$OUTPUT_GROUP" | "${process[@]}"
[[ -s "$PROCESSED" ]] || fail "GROMACS did not produce a non-empty $PROCESSED."

printf '\nStep 2: keep every %sth frame for VMD...\n' "$STRIDE"
printf '%s\n' "$OUTPUT_GROUP" | "${preview[@]}"
[[ -s "$PREVIEW" ]] || fail "GROMACS did not produce a non-empty $PREVIEW."
printf '\nFinished. Analyse %s; view %s in VMD.\n' "$PROCESSED" "$PREVIEW"
