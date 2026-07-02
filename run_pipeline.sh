#!/bin/bash
# =============================================================================
# run_pipeline.sh
# ONE COMMAND to run the complete ERY2-4/CTLA-4 docking pipeline
#
# Usage:
#   ./run_pipeline.sh                    # Full pipeline
#   ./run_pipeline.sh --step 1           # Run specific step
#   ./run_pipeline.sh --from 3           # Resume from step 3
#   ./run_pipeline.sh --dry-run          # Check setup only
#   ./run_pipeline.sh --method local     # Use local peptide PDB
#
# Reference: Ramanayake et al. ACS Chem. Biol. 2020, 15, 360-368
# PDB: 1I8L
# =============================================================================

set -euo pipefail

# ─── Colors ───────────────────────────────────────────────────────────────────
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# ─── Defaults ─────────────────────────────────────────────────────────────────
STEP=""
FROM_STEP=1
DRY_RUN=false
PEPTIDE_METHOD="pepfold3"
PEPTIDE_PDB=""
LOG_DIR="logs"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
MAIN_LOG="$LOG_DIR/pipeline_$TIMESTAMP.log"

# ─── Banner ───────────────────────────────────────────────────────────────────
print_banner() {
    echo -e "${BOLD}${BLUE}"
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║        ERY2-4 / CTLA-4 Automated Docking Pipeline           ║"
    echo "║                                                              ║"
    echo "║  Reference: Ramanayake et al. ACS Chem. Biol. 2020         ║"
    echo "║  PDB: 1I8L | Target: CTLA-4 | Peptide: ERY2-4 (HLH)       ║"
    echo "╚══════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
}

# ─── Logging ──────────────────────────────────────────────────────────────────
log() {
    local level=$1
    shift
    local msg="$*"
    local timestamp=$(date +"%H:%M:%S")
    case "$level" in
        INFO)  echo -e "${GREEN}[${timestamp}] ✓ ${msg}${NC}" | tee -a "$MAIN_LOG" ;;
        WARN)  echo -e "${YELLOW}[${timestamp}] ⚠ ${msg}${NC}" | tee -a "$MAIN_LOG" ;;
        ERROR) echo -e "${RED}[${timestamp}] ✗ ${msg}${NC}" | tee -a "$MAIN_LOG" ;;
        STEP)  echo -e "${CYAN}[${timestamp}] ► ${msg}${NC}" | tee -a "$MAIN_LOG" ;;
        *)     echo -e "[${timestamp}] ${msg}" | tee -a "$MAIN_LOG" ;;
    esac
}

# ─── Argument Parser ──────────────────────────────────────────────────────────
parse_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            --step)      STEP="$2";            shift 2 ;;
            --from)      FROM_STEP="$2";       shift 2 ;;
            --dry-run)   DRY_RUN=true;         shift   ;;
            --method)    PEPTIDE_METHOD="$2";  shift 2 ;;
            --peptide)   PEPTIDE_PDB="$2";     shift 2 ;;
            --help|-h)   show_help;            exit 0  ;;
            *)           log WARN "Unknown option: $1"; shift ;;
        esac
    done
}

show_help() {
    echo -e "${BOLD}Usage:${NC}"
    echo "  ./run_pipeline.sh [OPTIONS]"
    echo ""
    echo -e "${BOLD}Options:${NC}"
    echo "  --step N          Run only step N (1-6)"
    echo "  --from N          Resume pipeline from step N"
    echo "  --dry-run         Check dependencies only, don't run"
    echo "  --method METHOD   Peptide structure method: pepfold3|placeholder|local"
    echo "  --peptide PATH    Path to existing peptide PDB file"
    echo "  --help            Show this help"
    echo ""
    echo -e "${BOLD}Steps:${NC}"
    echo "  1 - Prepare CTLA-4 receptor (from PDB 1I8L)"
    echo "  2 - Prepare ERY2-4 peptide structure"
    echo "  3 - Generate HADDOCK3 restraints (AIR + disulfide)"
    echo "  4 - Run HADDOCK3 docking"
    echo "  5 - Analyze docking results"
    echo "  6 - Generate visualizations (PyMOL + plots)"
    echo ""
    echo -e "${BOLD}Examples:${NC}"
    echo "  ./run_pipeline.sh                          # Full pipeline"
    echo "  ./run_pipeline.sh --step 5                 # Analysis only"
    echo "  ./run_pipeline.sh --from 3                 # From restraints"
    echo "  ./run_pipeline.sh --peptide my_ery24.pdb   # Custom peptide"
    echo "  ./run_pipeline.sh --dry-run                # Check setup"
}

# ─── Dependency Checker ───────────────────────────────────────────────────────
check_dependencies() {
    log STEP "Checking dependencies..."
    local missing=0

    # Python
    if command -v python3 &>/dev/null; then
        PY_VER=$(python3 --version 2>&1)
        log INFO "Python: $PY_VER"
    else
        log ERROR "Python3 not found"
        ((missing++))
    fi

    # pdb-tools
    if command -v pdb_fetch &>/dev/null; then
        log INFO "pdb-tools: installed"
    else
        log WARN "pdb-tools: not found (will auto-install)"
    fi

    # HADDOCK3
    if command -v haddock3 &>/dev/null; then
        H3_VER=$(haddock3 --version 2>&1 || echo "unknown")
        log INFO "HADDOCK3: $H3_VER"
    else
        log WARN "HADDOCK3: not found"
        log WARN "  Install: cd haddock3 && pip install -e ."
        if [ "$STEP" = "4" ] || [ -z "$STEP" ]; then
            log ERROR "HADDOCK3 required for docking step"
            ((missing++))
        fi
    fi

    # Python packages
    python3 -c "import requests" 2>/dev/null && \
        log INFO "requests: installed" || \
        log WARN "requests: missing (pip install requests)"

    python3 -c "import numpy" 2>/dev/null && \
        log INFO "numpy: installed" || \
        log WARN "numpy: missing (pip install numpy)"

    python3 -c "import matplotlib" 2>/dev/null && \
        log INFO "matplotlib: installed" || \
        log WARN "matplotlib: missing (pip install matplotlib)"

    # PRODIGY
    if command -v prodigy &>/dev/null; then
        log INFO "PRODIGY: installed"
    else
        log WARN "PRODIGY: not found (pip install prodigy-prot)"
    fi

    # PyMOL
    if command -v pymol &>/dev/null; then
        log INFO "PyMOL: installed"
    else
        log WARN "PyMOL: not found (needed for visualization)"
    fi

    if [ "$missing" -gt 0 ]; then
        log ERROR "$missing critical dependency/dependencies missing"
        return 1
    fi

    log INFO "All critical dependencies satisfied"
    return 0
}

# ─── Environment Setup ────────────────────────────────────────────────────────
setup_environment() {
    log STEP "Setting up environment..."

    # Create directories
    mkdir -p logs \
             data/input \
             data/output/01_receptor \
             data/output/02_peptide \
             data/output/03_restraints \
             data/output/04_docking \
             data/output/05_analysis \
             data/output/06_visualization \
             results

    # Install missing Python packages
    pip install pdb-tools requests numpy matplotlib prodigy-prot -q \
        && log INFO "Python packages installed/updated" \
        || log WARN "Some packages may not have installed"

    log INFO "Environment ready"
}

# ─── Step Runner ──────────────────────────────────────────────────────────────
run_step() {
    local step_num=$1
    local step_name=$2
    local step_log="$LOG_DIR/step${step_num}_$TIMESTAMP.log"

    echo ""
    echo -e "${BOLD}${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    log STEP "STEP $step_num: $step_name"
    echo -e "${BOLD}${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

    local start_time=$(date +%s)

    case $step_num in
        1)
            bash scripts/01_prepare_receptor.sh 2>&1 | tee "$step_log"
            ;;
        2)
            local args="--method $PEPTIDE_METHOD"
            [ -n "$PEPTIDE_PDB" ] && args="$args --input-pdb $PEPTIDE_PDB"
            python3 scripts/02_prepare_peptide.py $args 2>&1 | tee "$step_log"
            ;;
        3)
            python3 scripts/03_generate_restraints.py 2>&1 | tee "$step_log"
            ;;
        4)
            bash scripts/04_run_docking.sh 2>&1 | tee "$step_log"
            ;;
        5)
            python3 scripts/05_analyze_results.py 2>&1 | tee "$step_log"
            ;;
        6)
            python3 scripts/06_visualize.py 2>&1 | tee "$step_log"
            ;;
        *)
            log ERROR "Unknown step: $step_num"
            return 1
            ;;
    esac

    local exit_code=${PIPESTATUS[0]}
    local end_time=$(date +%s)
    local elapsed=$((end_time - start_time))

    if [ $exit_code -eq 0 ]; then
        log INFO "Step $step_num completed in ${elapsed}s"
    else
        log ERROR "Step $step_num FAILED (exit code: $exit_code)"
        log ERROR "Check log: $step_log"
        return 1
    fi
}

# ─── Results Summary ──────────────────────────────────────────────────────────
print_summary() {
    echo ""
    echo -e "${BOLD}${GREEN}"
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║                  PIPELINE COMPLETE ✓                        ║"
    echo "╚══════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"

    echo -e "${BOLD}Output Files:${NC}"
    echo "  Receptor PDB    : data/output/01_receptor/ctla4_final.pdb"
    echo "  Peptide PDB     : data/output/02_peptide/ery2-4_final.pdb"
    echo "  Restraints      : data/output/03_restraints/"
    echo "  Docking models  : data/output/04_docking/"
    echo "  Analysis report : data/output/05_analysis/analysis_report.txt"
    echo "  PyMOL scripts   : data/output/06_visualization/*.pml"
    echo "  Contact heatmap : data/output/06_visualization/contact_heatmap.png"
    echo "  HADDOCK scores  : data/output/06_visualization/haddock_scores.png"
    echo ""
    echo -e "${BOLD}Next Steps:${NC}"
    echo "  1. Open PyMOL:  pymol data/output/06_visualization/01_overview.pml"
    echo "  2. View report: cat data/output/05_analysis/analysis_report.txt"
    echo "  3. Check logs:  ls $LOG_DIR/"
    echo ""
    echo -e "${BOLD}Experimental Reference (your paper):${NC}"
    echo "  KD  = 196.8 ± 2.3 nM  (SPR, 25°C)"
    echo "  IC50 = 1.1 ± 0.03 μM  (SPR inhibition)"
    echo "  IC50 = 631.9 nM        (Cell-based, DCs)"
    echo "  Trp33 (L33W) = primary anchor residue"
    echo ""
}

# ─── Main ─────────────────────────────────────────────────────────────────────
main() {
    parse_args "$@"

    mkdir -p "$LOG_DIR"
    print_banner

    log INFO "Pipeline started: $(date)"
    log INFO "Log file: $MAIN_LOG"

    # Dependency check
    check_dependencies || {
        log ERROR "Dependency check failed. Fix issues and retry."
        exit 1
    }

    if $DRY_RUN; then
        log INFO "Dry run complete. All checks passed."
        exit 0
    fi

    # Setup
    setup_environment

    # Determine which steps to run
    if [ -n "$STEP" ]; then
        # Run single step
        run_step "$STEP" "$(get_step_name $STEP)"
    else
        # Run all steps from FROM_STEP
        for step in $(seq "$FROM_STEP" 6); do
            run_step "$step" "$(get_step_name $step)" || {
                log ERROR "Pipeline failed at step $step"
                log ERROR "Resume with: ./run_pipeline.sh --from $step"
                exit 1
            }
        done
        print_summary
    fi
}

get_step_name() {
    case $1 in
        1) echo "Prepare CTLA-4 Receptor" ;;
        2) echo "Prepare ERY2-4 Peptide" ;;
        3) echo "Generate Restraints" ;;
        4) echo "Run HADDOCK3 Docking" ;;
        5) echo "Analyze Results" ;;
        6) echo "Generate Visualizations" ;;
        *) echo "Unknown" ;;
    esac
}

main "$@"
