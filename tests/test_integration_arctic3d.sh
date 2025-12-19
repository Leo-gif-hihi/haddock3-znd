#!/bin/bash

# Test script for a.sh with ARCTIC3D integration
# This script mocks arctic3d and arctic3d-restraints to verify the footprint matching logic.

set -e

# -----------------------------------------------------------------------------
# Setup Environment
# -----------------------------------------------------------------------------
TEST_DIR=$(mktemp -d)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
A_SH="$ROOT_DIR/script/a.sh"

echo "Running tests in $TEST_DIR"

# Create Mock Bin Directory
MOCK_BIN="$TEST_DIR/mock_bin"
mkdir -p "$MOCK_BIN"
export PATH="$MOCK_BIN:$PATH"

# -----------------------------------------------------------------------------
# Create Mocks
# -----------------------------------------------------------------------------

# Mock arctic3d
# Usage: arctic3d <uniprot_id> --pdb_to_use=<pdb_id> --chain=<chain>
cat <<'EOF' > "$MOCK_BIN/arctic3d"
#!/bin/bash
UNIPROT_ID="$1"
PDB_ID=""
CHAIN=""

for arg in "$@"; do
    if [[ "$arg" == --chain=* ]]; then
        CHAIN="${arg#*=}"
    elif [[ "$arg" == --pdb_to_use=* ]]; then
        PDB_ID="${arg#*=}"
    fi
done

echo "MOCK ARCTIC3D: Processing UniProt $UNIPROT_ID PDB $PDB_ID Chain $CHAIN"

# Create output directory structure expected by a.sh
# a.sh runs this inside chain_work_dir
# It expects a directory named "arctic3d-*" to be created
OUT_DIR="arctic3d-${UNIPROT_ID}"
mkdir -p "$OUT_DIR"

# Create a dummy PDB file inside the output directory
# This PDB should have the "Original" numbering (e.g., 101, 102, 103)
# Coordinates must match the input file (10.0, 20.0, 30.0)
cat <<PDB > "$OUT_DIR/${PDB_ID}.pdb"
ATOM      1  CA  ALA $CHAIN 101      10.000  10.000  10.000  1.00  0.00           C
ATOM      2  CA  VAL $CHAIN 102      20.000  20.000  20.000  1.00  0.00           C
ATOM      3  CA  GLY $CHAIN 103      30.000  30.000  30.000  1.00  0.00           C
PDB

echo "MOCK ARCTIC3D: Created $OUT_DIR/${PDB_ID}.pdb"
EOF
chmod +x "$MOCK_BIN/arctic3d"

# Mock arctic3d-restraints
# Usage: arctic3d-restraints --r1 ... --output <dir>
cat <<'EOF' > "$MOCK_BIN/arctic3d-restraints"
#!/bin/bash
OUTPUT_DIR=""
CH1="A"
CH2="B"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --ch1)
            CH1="$2"
            shift 2
            ;;
        --ch2)
            CH2="$2"
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done

echo "MOCK ARCTIC3D-RESTRAINTS: Output to $OUTPUT_DIR, Chains $CH1 $CH2"
mkdir -p "$OUTPUT_DIR"

# Create a dummy ambig.tbl with "Original" numbering (101, 103)
# We expect a.sh to translate this to (1, 3)
cat <<TBL > "$OUTPUT_DIR/mock_ambig.tbl"
! Mock Restraints
assign (segid $CH1 and resid 101) (segid $CH2 and resid 101) 2.0 2.0 0.0
assign (segid $CH1 and resid 103) (segid $CH2 and resid 103) 2.0 2.0 0.0
TBL
EOF
chmod +x "$MOCK_BIN/arctic3d-restraints"

# Mock haddock3 (just to pass the check)
cat <<'EOF' > "$MOCK_BIN/haddock3"
#!/bin/bash
echo "MOCK HADDOCK3"
EOF
chmod +x "$MOCK_BIN/haddock3"

# -----------------------------------------------------------------------------
# Create Input Data
# -----------------------------------------------------------------------------

DATA_DIR="$TEST_DIR/data"
mkdir -p "$DATA_DIR"

# Create "Original" PDBs with numbering 101-103
# 1A00.pdb (Partner A)
cat <<PDB > "$DATA_DIR/1A00.pdb"
HEADER    TEST PDB 1A00
DBREF  1A00 A    1   100  UNP    P12345   P12345_HUMAN     1    100
ATOM      1  CA  ALA A 101      10.000  10.000  10.000  1.00  0.00           C
ATOM      2  CA  VAL A 102      20.000  20.000  20.000  1.00  0.00           C
ATOM      3  CA  GLY A 103      30.000  30.000  30.000  1.00  0.00           C
END
PDB

# 1B00.pdb (Partner B) - Same coordinates/residues for simplicity
cat <<PDB > "$DATA_DIR/1B00.pdb"
HEADER    TEST PDB 1B00
DBREF  1B00 B    1   100  UNP    P67890   P67890_HUMAN     1    100
ATOM      1  CA  ALA B 101      10.000  10.000  10.000  1.00  0.00           C
ATOM      2  CA  VAL B 102      20.000  20.000  20.000  1.00  0.00           C
ATOM      3  CA  GLY B 103      30.000  30.000  30.000  1.00  0.00           C
END
PDB

# -----------------------------------------------------------------------------
# Run a.sh
# -----------------------------------------------------------------------------

echo "Running a.sh..."
# We use --dry-run to avoid running haddock3 itself, but a.sh should still generate restraints
if ! bash "$A_SH" \
    --partner A="$DATA_DIR/1A00.pdb" \
    --partner B="$DATA_DIR/1B00.pdb" \
    --use-arctic3d \
    --out "$TEST_DIR/result" \
    --project "test_run" \
    --dry-run; then
    
    echo "a.sh failed!"
    # Find and print any log files
    find "$TEST_DIR" -name "*.log" -print0 | xargs -0 -I {} sh -c 'echo "=== {} ==="; cat "{}"'
    exit 1
fi

# -----------------------------------------------------------------------------
# Verify Results
# -----------------------------------------------------------------------------

RESULT_DIR="$TEST_DIR/result/test_run"
PAIR_DIR="$RESULT_DIR/PAIR_A_vs_B"
RESTRAINTS_DIR="$PAIR_DIR/restraints"

echo "Checking results in $RESTRAINTS_DIR..."

# Check for ARCTIC3D restraints file
ARCTIC_TBL="$RESTRAINTS_DIR/arctic3d_PAIR_A_vs_B_ambig.tbl"
COMBINED_TBL="$RESTRAINTS_DIR/PAIR_A_vs_B_combined_ambig.tbl"

TARGET_TBL=""
if [[ -f "$COMBINED_TBL" ]]; then
    TARGET_TBL="$COMBINED_TBL"
elif [[ -f "$ARCTIC_TBL" ]]; then
    TARGET_TBL="$ARCTIC_TBL"
else
    echo "ERROR: No restraints file found (checked $COMBINED_TBL and $ARCTIC_TBL)"
    # Debug: list files
    ls -R "$RESULT_DIR"
    exit 1
fi

echo "Found restraints file: $TARGET_TBL. Content:"
cat "$TARGET_TBL"

# Check for translated residues
# Original was 101, 103.
# a.sh renumbers 101->1, 102->2, 103->3.
# So we expect "resid 1" and "resid 3".
# We do NOT expect "resid 101".

if grep -q "resid 1)" "$TARGET_TBL" && grep -q "resid 3)" "$TARGET_TBL"; then
    echo "SUCCESS: Found translated residues (1 and 3)."
else
    echo "FAILURE: Did not find expected translated residues."
    exit 1
fi

if grep -q "resid 101)" "$TARGET_TBL"; then
    echo "FAILURE: Found untranslated residue (101)."
    exit 1
fi

echo "Test Passed!"
rm -rf "$TEST_DIR"
