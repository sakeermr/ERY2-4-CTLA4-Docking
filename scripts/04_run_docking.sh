#!/bin/bash
# =============================================================================
# 04_run_docking.sh
# Run HADDOCK3 docking: ERY2-4 onto CTLA-4
# =============================================================================

set -euo pipefail

echo "=============================================="
echo " STEP 4: Running HADDOCK3 Docking"
echo " ERY2-4 (HLH peptide) -> CTLA-4 (1I8L)"
echo "=============================================="

CONFIG="config/haddock3.cfg"
OUTPUT_DIR="data/output/04_docking"
LOG_FILE="$OUTPUT_DIR/haddock3_run.log"

mkdir -p "$OUTPUT_DIR"

# --- Check HADDOCK3 installation ---
if ! command -v haddock3 &> /dev/null; then
    echo "[ERROR] HADDOCK3 not found."
    echo "        Install with: pip install haddock3"
    echo "        Or: cd haddock3 && pip install -e ."
    exit 1
fi

# --- Check input files ---
echo "[CHECK] Verifying input files..."
files=(
    "data/output/01_receptor/ctla4_final.pdb"
    "data/output/02_peptide/ery2-4_final.pdb"
    "data/output/03_restraints/air_restraints.tbl"
    "data/output/03_restraints/ssbond_restraints.tbl"
    "$CONFIG"
)

for f in "${files[@]}"; do
    if [ -f "$f" ]; then
        echo "      ✓ $f"
    else
        echo "      ✗ MISSING: $f"
        echo "[ERROR] Run previous steps first: make prepare"
        exit 1
    fi
done

echo ""
echo "[INFO] Starting HADDOCK3..."
echo "       Config  : $CONFIG"
echo "       Log     : $LOG_FILE"
echo "       Started : $(date)"
echo ""

# --- Run HADDOCK3 ---
haddock3 "$CONFIG" 2>&1 | tee "$LOG_FILE"

echo ""
echo "=============================================="
echo " HADDOCK3 Run Complete"
echo " Finished : $(date)"
echo "=============================================="
echo ""
echo " Results : $OUTPUT_DIR/ery2-4_ctla4_run/"
echo " Next    : make analyze"
