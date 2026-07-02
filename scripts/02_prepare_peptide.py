#!/usr/bin/env python3
"""
02_prepare_peptide.py
=====================
Prepare ERY2-4 peptide structure for HADDOCK3 docking.

Strategy:
  1. Submit sequence to PEP-FOLD3 API (online) OR
     use a locally installed structure predictor
  2. Apply disulfide bond constraint (Cys1 - Cys40)
  3. Rename chain to B
  4. Validate and save

ERY2-4 sequence: CAWGQAILEGELAWLEGGGGGAGQLADLKRQLAWWKQAC
Key features:
  - Disulfide: Cys1 -- Cys40 (cyclic HLH structure)
  - Critical: Trp33 (L33W mutation - primary anchor)
  - alpha-helical conformation (confirmed by CD spectroscopy)

Reference: Ramanayake et al. ACS Chem. Biol. 2020, 15, 360-368
"""

import os
import sys
import json
import time
import requests
import argparse
import subprocess
from pathlib import Path

# ─── Constants ────────────────────────────────────────────────────────────────

SEQUENCE     = "CAWGQAILEGELAWLEGGGGGAGQLADLKRQLAWWKQAC"
PEPTIDE_NAME = "ERY2-4"
SS_BOND      = (1, 40)   # Cys1 - Cys40 disulfide
OUTPUT_DIR   = Path("data/output/02_peptide")
INPUT_DIR    = Path("data/input")

PEPFOLD3_URL = "https://bioserv.rpbs.univ-paris-diderot.fr/services/PEP-FOLD3"

# ─── Helpers ──────────────────────────────────────────────────────────────────

def print_header():
    print("=" * 60)
    print(f"  STEP 2: Preparing ERY2-4 Peptide Structure")
    print("=" * 60)
    print(f"  Sequence : {SEQUENCE}")
    print(f"  Length   : {len(SEQUENCE)} aa")
    print(f"  Disulfide: Cys{SS_BOND[0]} -- Cys{SS_BOND[1]}")
    print(f"  Critical : Trp33 (L33W mutation)")
    print("=" * 60)


def submit_pepfold3(sequence: str) -> str:
    """
    Submit sequence to PEP-FOLD3 server.
    Returns job ID or raises on failure.
    """
    print("\n[INFO] Submitting to PEP-FOLD3...")
    payload = {
        "seq": sequence,
        "nb_run": 5,
        "name": PEPTIDE_NAME
    }
    try:
        r = requests.post(
            f"{PEPFOLD3_URL}/submit",
            json=payload,
            timeout=30
        )
        r.raise_for_status()
        job_id = r.json().get("job_id")
        print(f"[INFO] Job submitted. ID: {job_id}")
        return job_id
    except Exception as e:
        print(f"[WARN] PEP-FOLD3 submission failed: {e}")
        return None


def poll_pepfold3(job_id: str, max_wait: int = 300) -> bool:
    """Poll PEP-FOLD3 until job completes."""
    print(f"[INFO] Polling PEP-FOLD3 job {job_id}...")
    elapsed = 0
    while elapsed < max_wait:
        try:
            r = requests.get(
                f"{PEPFOLD3_URL}/status/{job_id}",
                timeout=15
            )
            status = r.json().get("status", "unknown")
            print(f"      Status: {status} ({elapsed}s elapsed)")
            if status == "completed":
                return True
            elif status == "failed":
                return False
        except Exception:
            pass
        time.sleep(15)
        elapsed += 15
    return False


def download_pepfold3_result(job_id: str, out_path: Path) -> bool:
    """Download best PEP-FOLD3 model."""
    try:
        r = requests.get(
            f"{PEPFOLD3_URL}/result/{job_id}/model_1.pdb",
            timeout=30
        )
        r.raise_for_status()
        out_path.write_text(r.text)
        print(f"[INFO] Downloaded model: {out_path}")
        return True
    except Exception as e:
        print(f"[WARN] Download failed: {e}")
        return False


def build_placeholder_pdb(sequence: str, out_path: Path):
    """
    Build a linear extended-chain PDB as fallback.
    Real docking should use PEP-FOLD3 / AlphaFold model.
    """
    print("[WARN] Building extended-chain placeholder PDB.")
    print("       Replace with PEP-FOLD3 / AlphaFold model for production!")

    lines = ["REMARK  ERY2-4 placeholder - replace with predicted structure\n"]
    aa3 = {
        'C':'CYS','A':'ALA','W':'TRP','G':'GLY','S':'SER',
        'I':'ILE','L':'LEU','E':'GLU','T':'THR','D':'ASP',
        'N':'ASN','K':'LYS','R':'ARG','Q':'GLN','M':'MET',
        'F':'PHE','Y':'TYR','P':'PRO','H':'HIS','V':'VAL'
    }

    atom_num = 1
    for i, aa in enumerate(sequence, start=1):
        res3 = aa3.get(aa, 'GLY')
        x = float(i) * 3.8
        y = 0.0
        z = 0.0
        line = (
            f"ATOM  {atom_num:5d}  CA  {res3} B{i:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
        )
        lines.append(line)
        atom_num += 1

    lines.append("END\n")
    out_path.write_text("".join(lines))
    print(f"[INFO] Placeholder PDB saved: {out_path}")


def add_ssbond_record(pdb_path: Path, res1: int, res2: int, chain: str = "B"):
    """Prepend SSBOND record to PDB file."""
    content = pdb_path.read_text()
    ssbond = (
        f"SSBOND   1 CYS {chain} {res1:4d}    CYS {chain} {res2:4d}"
        f"                       2.03\n"
    )
    pdb_path.write_text(ssbond + content)
    print(f"[INFO] SSBOND record added: CYS{res1} -- CYS{res2} (chain {chain})")


def rename_chain(pdb_path: Path, new_chain: str = "B"):
    """Rename all chains to B using pdb-tools."""
    try:
        result = subprocess.run(
            ["pdb_chain", f"-{new_chain}", str(pdb_path)],
            capture_output=True, text=True
        )
        if result.returncode == 0:
            pdb_path.write_text(result.stdout)
            print(f"[INFO] Chain renamed to: {new_chain}")
        else:
            print(f"[WARN] Chain rename failed: {result.stderr}")
    except FileNotFoundError:
        print("[WARN] pdb-tools not found. Chain rename skipped.")


def validate_pdb(pdb_path: Path):
    """Basic validation of PDB file."""
    content = pdb_path.read_text()
    atom_lines = [l for l in content.split('\n') if l.startswith('ATOM')]
    res_nums = sorted(set(int(l[22:26]) for l in atom_lines if l.strip()))

    print(f"\n[VALIDATION] {pdb_path.name}")
    print(f"  ATOM records  : {len(atom_lines)}")
    print(f"  Residues      : {len(res_nums)}")
    print(f"  Residue range : {min(res_nums)} - {max(res_nums)}")
    print(f"  Has SSBOND    : {'SSBOND' in content}")
    print(f"  Sequence len  : {len(SEQUENCE)}")

    # Check Cys positions
    cys_residues = [
        int(l[22:26]) for l in atom_lines
        if 'CYS' in l and ' CA ' in l
    ]
    print(f"  CYS positions : {cys_residues}")


def save_metadata(out_dir: Path, method: str):
    """Save preparation metadata as JSON."""
    meta = {
        "peptide": PEPTIDE_NAME,
        "sequence": SEQUENCE,
        "length": len(SEQUENCE),
        "disulfide": {"Cys1": SS_BOND[0], "Cys2": SS_BOND[1]},
        "structure_method": method,
        "chain": "B",
        "key_residues": {
            "Trp3":  "N-helix anchor",
            "Trp14": "loop region contact",
            "Trp33": "PRIMARY anchor - L33W critical mutation",
            "Trp39": "C-helix contact"
        },
        "reference": "Ramanayake et al. ACS Chem. Biol. 2020, 15, 360-368"
    }
    meta_path = out_dir / "peptide_metadata.json"
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"[INFO] Metadata saved: {meta_path}")


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Prepare ERY2-4 peptide for docking")
    parser.add_argument(
        "--method",
        choices=["pepfold3", "placeholder", "local"],
        default="pepfold3",
        help="Structure prediction method"
    )
    parser.add_argument(
        "--input-pdb",
        type=str,
        default=None,
        help="Path to existing peptide PDB (skip prediction)"
    )
    args = parser.parse_args()

    print_header()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    raw_pdb   = OUTPUT_DIR / "ery2-4_raw.pdb"
    final_pdb = OUTPUT_DIR / "ery2-4_final.pdb"
    method_used = args.method

    # ── Get structure ──────────────────────────────────────────────────────────
    if args.input_pdb:
        # User provides existing PDB
        import shutil
        shutil.copy(args.input_pdb, raw_pdb)
        print(f"[INFO] Using provided PDB: {args.input_pdb}")
        method_used = "user_provided"

    elif args.method == "pepfold3":
        job_id = submit_pepfold3(SEQUENCE)
        if job_id and poll_pepfold3(job_id):
            success = download_pepfold3_result(job_id, raw_pdb)
            if not success:
                build_placeholder_pdb(SEQUENCE, raw_pdb)
                method_used = "placeholder"
        else:
            print("[WARN] PEP-FOLD3 failed. Falling back to placeholder.")
            build_placeholder_pdb(SEQUENCE, raw_pdb)
            method_used = "placeholder"

    else:
        build_placeholder_pdb(SEQUENCE, raw_pdb)
        method_used = "placeholder"

    # ── Post-process ──────────────────────────────────────────────────────────
    import shutil
    shutil.copy(raw_pdb, final_pdb)

    rename_chain(final_pdb, new_chain="B")
    add_ssbond_record(final_pdb, SS_BOND[0], SS_BOND[1], chain="B")
    validate_pdb(final_pdb)
    save_metadata(OUTPUT_DIR, method_used)

    print(f"\n[DONE] ERY2-4 peptide prepared: {final_pdb}")
    print(f"       Method: {method_used}")
    if method_used == "placeholder":
        print("\n[!] WARNING: Using placeholder structure.")
        print("    For publication-quality results, use PEP-FOLD3 or AlphaFold.")
        print("    Submit manually at: https://bioserv.rpbs.univ-paris-diderot.fr/services/PEP-FOLD3/")


if __name__ == "__main__":
    main()
