# =============================================================================
# Makefile: ERY2-4 / CTLA-4 Docking Pipeline
# Reference: Ramanayake et al. ACS Chem. Biol. 2020, 15, 360-368
# =============================================================================

.PHONY: all prepare receptor peptide restraints dock analyze visualize \
        clean clean-all install check help

# ── Colors ────────────────────────────────────────────────────────────────────
CYAN  = \033[0;36m
GREEN = \033[0;32m
BOLD  = \033[1m
NC    = \033[0m

# ── Default: full pipeline ────────────────────────────────────────────────────
all: check install prepare dock analyze visualize
	@echo -e "$(GREEN)$(BOLD)Pipeline complete!$(NC)"

# ── Installation ──────────────────────────────────────────────────────────────
install:
	@echo -e "$(CYAN)Installing dependencies...$(NC)"
	pip install pdb-tools requests numpy matplotlib prodigy-prot -q
	@echo -e "$(GREEN)Dependencies installed$(NC)"

# ── Dependency check ──────────────────────────────────────────────────────────
check:
	@echo -e "$(CYAN)Checking dependencies...$(NC)"
	@./run_pipeline.sh --dry-run

# ── Structure preparation ─────────────────────────────────────────────────────
prepare: receptor peptide restraints

receptor:
	@echo -e "$(CYAN)Step 1: Preparing CTLA-4 receptor...$(NC)"
	@bash scripts/01_prepare_receptor.sh

peptide:
	@echo -e "$(CYAN)Step 2: Preparing ERY2-4 peptide...$(NC)"
	@python3 scripts/02_prepare_peptide.py --method pepfold3

peptide-local:
	@echo -e "$(CYAN)Step 2: Building placeholder peptide...$(NC)"
	@python3 scripts/02_prepare_peptide.py --method placeholder

restraints:
	@echo -e "$(CYAN)Step 3: Generating restraints...$(NC)"
	@python3 scripts/03_generate_restraints.py

# ── Docking ───────────────────────────────────────────────────────────────────
dock:
	@echo -e "$(CYAN)Step 4: Running HADDOCK3...$(NC)"
	@bash scripts/04_run_docking.sh

# ── Analysis ──────────────────────────────────────────────────────────────────
analyze:
	@echo -e "$(CYAN)Step 5: Analyzing results...$(NC)"
	@python3 scripts/05_analyze_results.py

# ── Visualization ─────────────────────────────────────────────────────────────
visualize:
	@echo -e "$(CYAN)Step 6: Generating visualizations...$(NC)"
	@python3 scripts/06_visualize.py

# ── PyMOL launchers ───────────────────────────────────────────────────────────
view-overview:
	pymol data/output/06_visualization/01_overview.pml

view-trp33:
	pymol data/output/06_visualization/02_trp33_closeup.pml

view-comparison:
	pymol data/output/06_visualization/03_comparison.pml

# ── Resume from specific step ─────────────────────────────────────────────────
from-restraints:
	@./run_pipeline.sh --from 3

from-analysis:
	@./run_pipeline.sh --from 5

# ── Report ────────────────────────────────────────────────────────────────────
report:
	@cat data/output/05_analysis/analysis_report.txt

# ── Clean ─────────────────────────────────────────────────────────────────────
clean:
	@echo "Cleaning output files..."
	@rm -rf data/output/
	@rm -rf logs/
	@echo "Done."

clean-all: clean
	@echo "Cleaning all generated files..."
	@rm -rf data/input/1I8L.pdb
	@echo "Done."

# ── Help ──────────────────────────────────────────────────────────────────────
help:
	@echo ""
	@echo -e "$(BOLD)ERY2-4 / CTLA-4 Docking Pipeline$(NC)"
	@echo -e "Reference: Ramanayake et al. ACS Chem. Biol. 2020"
	@echo ""
	@echo -e "$(BOLD)Quick Start:$(NC)"
	@echo "  make all              # Full pipeline"
	@echo "  make check            # Check dependencies"
	@echo "  make install          # Install packages"
	@echo ""
	@echo -e "$(BOLD)Individual Steps:$(NC)"
	@echo "  make receptor         # Step 1: Prepare CTLA-4"
	@echo "  make peptide          # Step 2: Prepare ERY2-4"
	@echo "  make restraints       # Step 3: Generate restraints"
	@echo "  make dock             # Step 4: Run HADDOCK3"
	@echo "  make analyze          # Step 5: Analyze results"
	@echo "  make visualize        # Step 6: Generate figures"
	@echo ""
	@echo -e "$(BOLD)Visualization:$(NC)"
	@echo "  make view-overview    # Open PyMOL overview"
	@echo "  make view-trp33       # Open Trp33 close-up"
	@echo "  make view-comparison  # Open B7-1 vs ERY2-4"
	@echo "  make report           # Print analysis report"
	@echo ""
	@echo -e "$(BOLD)Utilities:$(NC)"
	@echo "  make clean            # Remove output files"
	@echo "  make clean-all        # Remove everything"
