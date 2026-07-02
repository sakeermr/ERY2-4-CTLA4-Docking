#!/usr/bin/env python3
"""
06_visualize.py
===============
Generate PyMOL visualization scripts and publication-quality figures
for ERY2-4 / CTLA-4 docking results.

Outputs:
  1. PyMOL .pml session script
  2. Contact heatmap (matplotlib)
  3. Binding site comparison figure (ERY2-4 vs B7-1)
  4. Trp33 close-up view script

Reference: Ramanayake et al. ACS Chem. Biol. 2020, 15, 360-368
"""

import json
import numpy as np
from pathlib import Path

# ─── Paths ────────────────────────────────────────────────────────────────────
ANALYSIS_DIR  = Path("data/output/05_analysis")
DOCKING_DIR   = Path("data/output/04_docking/ery2-4_ctla4_run")
OUTPUT_DIR    = Path("data/output/06_visualization")
REFERENCE_PDB = Path("data/input/1I8L.pdb")

# ─── Color Scheme ─────────────────────────────────────────────────────────────
COLORS = {
    "ctla4":          "slate",
    "b71":            "wheat",
    "ery24":          "tv_green",
    "mypppy":         "red",
    "trp33":          "yellow",
    "interface":      "orange",
    "b71_shared":     "magenta",
    "ery24_unique":   "cyan",
    "background":     "white"
}

# ─── PyMOL Script Generator ───────────────────────────────────────────────────

def generate_pymol_overview(best_model: str = None) -> str:
    """
    Generate PyMOL script for overview of ERY2-4/CTLA-4 complex.
    """
    model_path = best_model or "data/output/04_docking/best_model.pdb"

    script = f"""# ============================================================
# PyMOL Visualization: ERY2-4 / CTLA-4 Docking
# Reference: Ramanayake et al. ACS Chem. Biol. 2020
# PDB: 1I8L
# ============================================================

# --- Setup ---
bg_color {COLORS['background']}
set ray_shadows, 0
set ray_opaque_background, off
set antialias, 2
set hash_max, 300

# --- Load structures ---
load {REFERENCE_PDB}, reference_1I8L
load {model_path}, ery24_ctla4_complex

# --- Reference structure (1I8L) ---
select ctla4_ref,  reference_1I8L and chain C
select b71_ref,    reference_1I8L and chain A

color {COLORS['ctla4']},  ctla4_ref
color {COLORS['b71']},    b71_ref

show cartoon, ctla4_ref
show cartoon, b71_ref
set cartoon_transparency, 0.5, b71_ref

# --- Docked complex ---
select ctla4_docked,  ery24_ctla4_complex and chain C
select ery24_docked,  ery24_ctla4_complex and chain B

color {COLORS['ctla4']},   ctla4_docked
color {COLORS['ery24']},   ery24_docked

show cartoon, ctla4_docked
show cartoon, ery24_docked

# --- Align CTLA-4 chains ---
align ctla4_docked, ctla4_ref

# --- MYPPPY motif (97-102) ---
select mypppy, (reference_1I8L or ery24_ctla4_complex) and chain C and resi 97-102
color {COLORS['mypppy']}, mypppy
show sticks, mypppy
label mypppy and name CA, "MYPPPY(%s%s)" % (resn, resi)

# --- Key CTLA-4 interface residues ---
select ctla4_interface, chain C and resi 38+48+50+93+97+98+99+100+101+102+110
color {COLORS['interface']}, ctla4_interface and ery24_ctla4_complex
show sticks, ctla4_interface and ery24_ctla4_complex

# --- ERY2-4 key residues ---
select ery24_trp,   ery24_docked and resn TRP
select ery24_trp33, ery24_docked and resi 33
select ery24_cys,   ery24_docked and resn CYS

color {COLORS['trp33']}, ery24_trp33
show sticks, ery24_trp
show spheres, ery24_trp33

# --- Disulfide bond ---
select ss_bond, ery24_docked and resi 1+40 and resn CYS
show sticks, ss_bond
color orange, ss_bond
distance ss_distance, ery24_docked and resi 1 and name SG, ery24_docked and resi 40 and name SG

# --- Labels ---
label ery24_trp33 and name CA, "Trp33 (L33W)"
label ctla4_docked and resi 97 and name CA, "Met97"
label ctla4_docked and resi 48 and name CA, "Tyr48"
label ctla4_docked and resi 102 and name CA, "Tyr102"

# --- Surface representation ---
create ctla4_surface, ctla4_docked
show surface, ctla4_surface
color {COLORS['ctla4']}, ctla4_surface
set transparency, 0.4, ctla4_surface

# --- B7-1 footprint on CTLA-4 ---
select b71_footprint, ctla4_docked and resi 38+48+50+93+97+98+99+100+101+102+110
color {COLORS['b71_shared']}, b71_footprint and ctla4_surface

# --- Contacts / H-bonds ---
distance hbonds, ery24_docked, ctla4_docked, 3.5, mode=2
color green, hbonds
hide labels, hbonds

# --- View settings ---
set_view (\\
     0.68,  0.45, -0.58,\\
    -0.17,  0.87,  0.46,\\
     0.71, -0.20,  0.67,\\
     0.00,  0.00, -80.0,\\
     0.00,  0.00,   0.00,\\
    50.00, 120.00, -20.00)

# --- Zoom ---
zoom ery24_ctla4_complex, 5
center ery24_trp33

# --- Save session ---
save data/output/06_visualization/ery24_ctla4_overview.pse

print "Overview scene saved."
"""
    return script


def generate_pymol_trp33_closeup() -> str:
    """
    PyMOL script for Trp33 close-up view.
    Shows Trp33 insertion into CTLA-4 hydrophobic pocket.
    """
    script = f"""# ============================================================
# PyMOL Close-up: Trp33 in CTLA-4 Hydrophobic Pocket
# ============================================================

reinitialize

bg_color {COLORS['background']}
set ray_shadows, 0

load data/output/04_docking/ery2-4_ctla4_run/best_model.pdb, complex

# --- Selections ---
select trp33,          complex and chain B and resi 33
select ctla4_pocket,   complex and chain C and resi 48+93+97+98
select mypppy,         complex and chain C and resi 97+98+99+100+101+102
select ery24_helix2,   complex and chain B and resi 28-40

# --- Colors ---
color {COLORS['ery24']},     complex and chain B
color {COLORS['ctla4']},     complex and chain C
color {COLORS['trp33']},     trp33
color {COLORS['mypppy']},    mypppy
color salmon,             ctla4_pocket

# --- Show sticks for key residues ---
show sticks, trp33
show sticks, ctla4_pocket
show sticks, mypppy
show cartoon, ery24_helix2
show cartoon, complex and chain C

# --- Surface of CTLA-4 pocket ---
create pocket_surface, ctla4_pocket
show surface, pocket_surface
set transparency, 0.3, pocket_surface
color salmon, pocket_surface

# --- Measure distances ---
distance trp33_met97, (complex and chain B and resi 33 and name CZ2), \\
                      (complex and chain C and resi 97 and name SD)

distance trp33_tyr48, (complex and chain B and resi 33 and name CE2), \\
                      (complex and chain C and resi 48 and name OH)

# --- Labels ---
label trp33 and name CZ2,        "Trp33 (L33W)"
label ctla4_pocket and resi 97 and name SD,  "Met97"
label ctla4_pocket and resi 48 and name OH,  "Tyr48"
label mypppy and resi 99 and name CA,        "Pro99"

# --- Zoom into Trp33 ---
zoom trp33, 8
set_view (\\
     0.90, -0.20,  0.39,\\
     0.10,  0.97,  0.23,\\
    -0.42, -0.15,  0.90,\\
     0.00,  0.00, -25.0,\\
     0.00,  0.00,   0.00,\\
    15.00,  40.00, -20.00)

save data/output/06_visualization/trp33_closeup.pse
print "Trp33 close-up saved."
"""
    return script


def generate_pymol_comparison() -> str:
    """
    PyMOL script comparing ERY2-4 vs B7-1 binding on CTLA-4.
    Side-by-side view.
    """
    script = f"""# ============================================================
# PyMOL Comparison: ERY2-4 vs B7-1 on CTLA-4
# Shows competitive binding to same epitope
# ============================================================

reinitialize
bg_color {COLORS['background']}

# --- Load reference (B7-1/CTLA-4) ---
load data/input/1I8L.pdb, b71_ctla4

# --- Load docked complex ---
load data/output/04_docking/ery2-4_ctla4_run/best_model.pdb, ery24_ctla4

# --- Align on CTLA-4 ---
align ery24_ctla4 and chain C, b71_ctla4 and chain C

# --- B7-1 complex ---
select b71,       b71_ctla4 and chain A
select ctla4_b71, b71_ctla4 and chain C

color {COLORS['b71']},    b71
color {COLORS['ctla4']},  ctla4_b71
show cartoon, b71
show cartoon, ctla4_b71

# --- ERY2-4 complex ---
select ery24,       ery24_ctla4 and chain B
select ctla4_ery24, ery24_ctla4 and chain C

color {COLORS['ery24']},  ery24
color {COLORS['ctla4']},  ctla4_ery24
show cartoon, ery24
show cartoon, ctla4_ery24

# --- MYPPPY on both ---
select mypppy_b71,   b71_ctla4 and chain C and resi 97-102
select mypppy_ery24, ery24_ctla4 and chain C and resi 97-102

color {COLORS['mypppy']}, mypppy_b71
color {COLORS['mypppy']}, mypppy_ery24
show sticks, mypppy_b71
show sticks, mypppy_ery24

# --- Trp33 of ERY2-4 ---
select trp33, ery24_ctla4 and chain B and resi 33
color {COLORS['trp33']}, trp33
show sticks, trp33

# --- Surfaces ---
create surface_b71,   ctla4_b71
create surface_ery24, ctla4_ery24
show surface, surface_b71
show surface, surface_ery24
set transparency, 0.5

color {COLORS['b71_shared']}, surface_b71 and resi 38+48+50+93+97+98+99+100+101+102+110
color {COLORS['b71_shared']}, surface_ery24 and resi 38+48+50+93+97+98+99+100+101+102+110

# --- Title labels ---
pseudoatom label_b71,   pos=[-20, 0, 0]
pseudoatom label_ery24, pos=[ 20, 0, 0]
label label_b71,   "B7-1/CTLA-4"
label label_ery24, "ERY2-4/CTLA-4"

zoom all, 10
save data/output/06_visualization/comparison_b71_vs_ery24.pse
print "Comparison scene saved."
"""
    return script


def generate_contact_heatmap(analysis_results: dict, out_dir: Path):
    """
    Generate a contact heatmap using matplotlib.
    CTLA-4 residues (x) vs ERY2-4 residues (y).
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors

        contact_pairs = analysis_results.get("interface", {}).get("contact_pairs", [])

        if not contact_pairs:
            print("[WARN] No contact pairs found for heatmap.")
            return

        # Build contact matrix
        rec_res = sorted(set(p["receptor_res"] for p in contact_pairs))
        lig_res = sorted(set(p["ligand_res"]   for p in contact_pairs))

        matrix = np.zeros((len(lig_res), len(rec_res)))
        lig_idx = {r: i for i, r in enumerate(lig_res)}
        rec_idx = {r: i for i, r in enumerate(rec_res)}

        for p in contact_pairs:
            r = rec_idx[p["receptor_res"]]
            l = lig_idx[p["ligand_res"]]
            # Use minimum distance
            if matrix[l, r] == 0 or p["distance"] < matrix[l, r]:
                matrix[l, r] = p["distance"]

        # Replace zeros with NaN
        matrix[matrix == 0] = np.nan

        # Plot
        fig, ax = plt.subplots(figsize=(14, 8))

        cmap = plt.cm.RdYlGn_r
        cmap.set_bad('white')

        im = ax.imshow(matrix, cmap=cmap, aspect='auto',
                      vmin=2.0, vmax=CONTACT_CUTOFF)

        # Axes
        ax.set_xticks(range(len(rec_res)))
        ax.set_xticklabels([f"{r}" for r in rec_res], rotation=90, fontsize=8)
        ax.set_yticks(range(len(lig_res)))
        ax.set_yticklabels([f"{r}" for r in lig_res], fontsize=8)

        ax.set_xlabel("CTLA-4 Residue Number", fontsize=12, fontweight='bold')
        ax.set_ylabel("ERY2-4 Residue Number", fontsize=12, fontweight='bold')
        ax.set_title(
            "ERY2-4 / CTLA-4 Contact Map\n"
            "Color = Distance (Å) | Green = Close contact | Red = Distant",
            fontsize=13, fontweight='bold'
        )

        # Highlight MYPPPY
        mypppy = [97, 98, 99, 100, 101, 102]
        for res in mypppy:
            if res in rec_idx:
                ax.axvline(x=rec_idx[res], color='red', alpha=0.3, linewidth=2)

        # Highlight Trp33
        if 33 in lig_idx:
            ax.axhline(y=lig_idx[33], color='gold', alpha=0.5, linewidth=3,
                      label='Trp33 (L33W)')

        # Colorbar
        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label("Distance (Å)", fontsize=11)

        # Annotations
        ax.text(0.02, 0.98,
                "Red lines = MYPPPY motif (97-102)\nGold line = Trp33 (L33W critical)",
                transform=ax.transAxes, fontsize=9,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        plt.tight_layout()
        heatmap_path = out_dir / "contact_heatmap.png"
        plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"[INFO] Contact heatmap saved: {heatmap_path}")

    except ImportError:
        print("[WARN] matplotlib not available. Install: pip install matplotlib")
    except Exception as e:
        print(f"[WARN] Heatmap generation failed: {e}")


def generate_score_plot(analysis_results: dict, out_dir: Path):
    """Plot HADDOCK scores distribution."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        scores = analysis_results.get("haddock_scores", [])
        if not scores:
            print("[WARN] No HADDOCK scores for plotting.")
            return

        score_vals = []
        for s in scores:
            try:
                score_vals.append(float(s.get("score", 0)))
            except (ValueError, TypeError):
                pass

        if not score_vals:
            return

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        # Score distribution
        ax1.hist(score_vals, bins=20, color='steelblue', edgecolor='black', alpha=0.7)
        ax1.axvline(x=min(score_vals), color='red', linewidth=2, label=f'Best: {min(score_vals):.1f}')
        ax1.set_xlabel("HADDOCK Score (a.u.)", fontsize=12)
        ax1.set_ylabel("Count", fontsize=12)
        ax1.set_title("HADDOCK Score Distribution", fontsize=13, fontweight='bold')
        ax1.legend()

        # Score vs model rank
        ax2.plot(range(1, len(score_vals)+1), sorted(score_vals),
                'o-', color='steelblue', markersize=4)
        ax2.set_xlabel("Model Rank", fontsize=12)
        ax2.set_ylabel("HADDOCK Score (a.u.)", fontsize=12)
        ax2.set_title("Score vs Rank", fontsize=13, fontweight='bold')
        ax2.axhline(y=sorted(score_vals)[0], color='red', linestyle='--',
                   label='Best model')
        ax2.legend()

        plt.suptitle("ERY2-4 / CTLA-4 Docking - HADDOCK Scores",
                    fontsize=14, fontweight='bold')
        plt.tight_layout()

        score_path = out_dir / "haddock_scores.png"
        plt.savefig(score_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"[INFO] Score plot saved: {score_path}")

    except ImportError:
        print("[WARN] matplotlib not available.")
    except Exception as e:
        print(f"[WARN] Score plot failed: {e}")


def main():
    print("=" * 60)
    print("  STEP 6: Generating Visualizations")
    print("=" * 60)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ── Load analysis results ──────────────────────────────────────────────────
    results_file = ANALYSIS_DIR / "analysis_results.json"
    analysis_results = {}
    if results_file.exists():
        with open(results_file) as f:
            analysis_results = json.load(f)
        print(f"[INFO] Loaded analysis results: {results_file}")
    else:
        print(f"[WARN] No analysis results found at {results_file}")
        print("       Run step 5 first: python scripts/05_analyze_results.py")

    # ── Find best docked model ─────────────────────────────────────────────────
    model_files = list(DOCKING_DIR.glob("**/*.pdb")) if DOCKING_DIR.exists() else []
    best_model = str(model_files[0]) if model_files else str(REFERENCE_PDB)

    # ── Generate PyMOL scripts ─────────────────────────────────────────────────
    scripts = {
        "01_overview.pml":    generate_pymol_overview(best_model),
        "02_trp33_closeup.pml": generate_pymol_trp33_closeup(),
        "03_comparison.pml":  generate_pymol_comparison()
    }

    for fname, content in scripts.items():
        path = OUTPUT_DIR / fname
        path.write_text(content)
        print(f"[INFO] PyMOL script: {path}")

    # ── Generate figures ───────────────────────────────────────────────────────
    print("\n[INFO] Generating contact heatmap...")
    generate_contact_heatmap(analysis_results, OUTPUT_DIR)

    print("[INFO] Generating score plot...")
    generate_score_plot(analysis_results, OUTPUT_DIR)

    # ── Usage instructions ─────────────────────────────────────────────────────
    instructions = """
HOW TO USE PyMOL SCRIPTS
=========================

1. Overview scene:
   pymol data/output/06_visualization/01_overview.pml

2. Trp33 close-up:
   pymol data/output/06_visualization/02_trp33_closeup.pml

3. ERY2-4 vs B7-1 comparison:
   pymol data/output/06_visualization/03_comparison.pml

Key visualization points:
  - CTLA-4  = slate blue cartoon
  - ERY2-4  = green cartoon
  - B7-1    = wheat cartoon
  - MYPPPY  = red sticks (97-102)
  - Trp33   = yellow spheres (L33W primary anchor)
  - Shared B7-1/ERY2-4 surface = magenta

For publication figures:
  PyMOL > ray 2400, 1800
  PyMOL > png figure.png, dpi=300
"""
    instr_path = OUTPUT_DIR / "HOW_TO_VISUALIZE.txt"
    instr_path.write_text(instructions)
    print(instructions)

    print(f"[DONE] Visualization files saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
