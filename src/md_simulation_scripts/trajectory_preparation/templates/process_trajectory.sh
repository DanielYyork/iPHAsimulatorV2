#!/usr/bin/env bash
# Read settings from instrcution.txt beside this script; requires only Bash and GROMACS.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

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
while IFS= read -r line || [[ -n "$line" ]]; do
    line="${line%$'\r'}"
    [[ "$line" =~ ^[[:space:]]*$ || "$line" =~ ^[[:space:]]*# ]] && continue
    [[ "$line" == *=* ]] || fail "Expected KEY=value in instrcution.txt: $line"
    key="${line%%=*}"
    value="${line#*=}"
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
for file in "$INPUT" "$TPR"; do
    [[ -r "$file" && -f "$file" ]] || fail "Missing or unreadable input: $file"
done
[[ "$TPR" == *.tpr ]] || fail "Use a matching production .tpr for -pbc mol."
[[ "$PROCESSED" != "$PREVIEW" ]] || fail "The two output filenames must differ."
for file in "$PROCESSED" "$PREVIEW"; do
    [[ "$file" == *.xtc && "$file" != */* ]] || fail "Outputs must be .xtc filenames inside this folder."
    [[ ! -e "$file" && ! -L "$file" ]] || fail "Output already exists: $file. Rename it or change the output setting."
done

# Generate a separate default index from the matching production TPR.
# Never modify the original simulation index.
[[ "$INDEX" == *.ndx && "$INDEX" != */* ]] || fail "INDEX must be a local .ndx filename, e.g. protein_center.ndx."
[[ ! -L "$INDEX" ]] || fail "INDEX must not be a symbolic link."
check_group() {
    awk -v wanted="$1" '
        /^[[:space:]]*\[/ {
            header=$0
            sub(/^[[:space:]]*\[[[:space:]]*/, "", header)
            sub(/[[:space:]]*\][[:space:]]*$/, "", header)
            selected=(header==wanted)
            next
        }
        selected && /^[[:space:]]*[0-9]/ { count+=NF }
        END { exit(count>0 ? 0 : 1) }
    ' "$INDEX" || fail "Missing or empty group $1 in $INDEX."
}
if [[ -e "$INDEX" ]]; then
    [[ -r "$INDEX" && -f "$INDEX" ]] || fail "Cannot read $INDEX."
    check_group "$CENTER_GROUP"
    check_group "$OUTPUT_GROUP"
fi
make_index=("$GMX" make_ndx -f "$TPR" -o "$INDEX")

process=("$GMX" trjconv -f "$INPUT" -s "$TPR" -n "$INDEX"
         -pbc mol -ur compact -center -o "$PROCESSED")
preview=("$GMX" trjconv -f "$PROCESSED" -s "$TPR" -n "$INDEX"
         -skip "$STRIDE" -o "$PREVIEW")

printf 'Full-resolution analysis trajectory: %s\nVMD preview: %s\n' "$PROCESSED" "$PREVIEW"
if "$DRY_RUN"; then
    if [[ ! -e "$INDEX" ]]; then
        printf '\nCreate index first (selection: q):\n'
        printf '%q ' "${make_index[@]}"; printf '\n'
    fi
    printf '\nStep 1 selections: %s, then %s\n' "$CENTER_GROUP" "$OUTPUT_GROUP"
    printf '%q ' "${process[@]}"; printf '\n'
    printf '\nStep 2 selection: %s\n' "$OUTPUT_GROUP"
    printf '%q ' "${preview[@]}"; printf '\n'
    printf '\nDry run only: no GROMACS commands executed or outputs written.\n'
    exit 0
fi

if [[ ! -e "$INDEX" ]]; then
    printf '\nCreating %s with default Protein and System groups...\n' "$INDEX"
    printf 'q\n' | "${make_index[@]}"
    check_group "$CENTER_GROUP"
    check_group "$OUTPUT_GROUP"
fi

printf '\nStep 1: centre on %s and compact-wrap all frames...\n' "$CENTER_GROUP"
printf '%s\n%s\n' "$CENTER_GROUP" "$OUTPUT_GROUP" | "${process[@]}"
[[ -s "$PROCESSED" ]] || fail "GROMACS did not produce a non-empty $PROCESSED."

printf '\nStep 2: keep every %sth frame for VMD...\n' "$STRIDE"
printf '%s\n' "$OUTPUT_GROUP" | "${preview[@]}"
[[ -s "$PREVIEW" ]] || fail "GROMACS did not produce a non-empty $PREVIEW."
printf '\nFinished. Analyse %s; view %s in VMD.\n' "$PROCESSED" "$PREVIEW"
