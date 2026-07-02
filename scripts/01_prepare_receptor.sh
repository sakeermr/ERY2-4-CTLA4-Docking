#!/bin/bash
# =============================================================================
# 01_prepare_receptor.sh
# Prepare CTLA-4 receptor from PDB 1I8L
# Target: CTLA-4 (chains C/D) from Human B7-1/CTLA-4 complex
# Reference: Stamper et al. Nature 2001, PDB: 1I8L
# =============================================================================

set -euo pipefail

echo "=============================================="
echo " STEP 1: Preparing CTLA-4 Receptor (1I8L)"
echo "=============================================="

# --- Paths ---
INPUT_DIR="data/input"
OUTPUT_DIR="data/output/01_receptor"
mkdir -p "$OUTPUT_DIR"

# --- Check pdb-tools ---
if ! command -v pdb_fetch &> /dev/null; then
    echo "[INFO] Installing pdb-tools..."
    pip install pdb-tools -q
fi

# --- Download 1I8L ---
echo "[1/6] Downloading PDB 1I8L..."
if [ ! -f "$INPUT_DIR/1I8L.pdb" ]; then
    pdb_fetch 1I8L > "$INPUT_DIR/1I8L.pdb"
    echo "      Downloaded: $INPUT_DIR/1I8L.pdb"
else
    echo "      Already exists: $INPUT_DIR/1I8L.pdb"
fi

# --- Extract CTLA-4 chain C ---
echo "[2/6] Extracting CTLA-4 chain C..."
pdb_selchain -C "$INPUT_DIR/1I8L.pdb" > "$OUTPUT_DIR/ctla4_chainC_raw.pdb"

# --- Remove HETATM (NAG, MAN glycans) ---
echo "[3/6] Removing HETATM records (glycans)..."
pdb_delhetatm "$OUTPUT_DIR/ctla4_chainC_raw.pdb" > "$OUTPUT_DIR/ctla4_nohet.pdb"

# --- Remove water molecules ---
echo "[4/6] Removing water molecules..."
grep -v "^HETATM\|HOH\|WAT" "$OUTPUT_DIR/ctla4_nohet.pdb" > "$OUTPUT_DIR/ctla4_nowat.pdb"

# --- Fix residue numbering ---
echo "[5/6] Fixing residue numbering..."
pdb_reres "$OUTPUT_DIR/ctla4_nowat.pdb" > "$OUTPUT_DIR/ctla4_renumbered.pdb"

# --- Tidy and validate ---
echo "[6/6] Final tidy and validation..."
pdb_tidy "$OUTPUT_DIR/ctla4_renumbered.pdb" > "$OUTPUT_DIR/ctla4_final.pdb"

echo ""
echo "[DONE] CTLA-4 receptor prepared: $OUTPUT_DIR/ctla4_final.pdb"
echo ""

# --- Quick stats ---
RESIDUE_COUNT=$(grep "^ATOM" "$OUTPUT_DIR/ctla4_final.pdb" | awk '{print $6}' | sort -u | wc -l)
echo "      Residue count: $RESIDUE_COUNT"
echo "      Chain: C"
echo "      Key binding residues: 38,48,50,93,97-102,110 (MYPPPY motif)"
