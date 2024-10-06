# Variables
SLURM_FILES := $(filter-out jobs/plot.slurm, $(wildcard jobs/*.slurm))
PLOT_SLURM := jobs/plot.slurm

# Default target: execute all SLURM job files except plot.slurm
all: execute

# Target to execute all SLURM job files except plot.slurm
execute: $(SLURM_FILES)
	@echo "Executing SLURM job scripts..."
	@for job in $(SLURM_FILES); do \
		echo "Submitting $$job..."; \
		sbatch $$job; \
	done

# Target to execute plot.slurm specifically
plot: $(PLOT_SLURM)
	@echo "Submitting plot SLURM job..."
	sbatch $(PLOT_SLURM)

# Target to clean any generated files (you can specify specific files if needed)
clean:
	@echo "Cleaning generated files..."
	rm -rf csv logs plots

# Phony targets
.PHONY: all execute plot clean
