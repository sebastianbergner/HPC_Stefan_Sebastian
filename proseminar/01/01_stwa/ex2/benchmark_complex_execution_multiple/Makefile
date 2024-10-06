# Variables
SLURM_FILES := $(sort $(filter-out jobs/plot.slurm, $(wildcard jobs/*.slurm)))
PLOT_SLURM := jobs/plot.slurm
REPS ?= 1  # Default repetitions if not specified

# Default target: execute all SLURM job files except plot.slurm
all: execute

# Target to clean and then execute all SLURM job files except plot.slurm, with repetitions
execute: clean $(SLURM_FILES)
	@echo "Executing SLURM job scripts $(REPS) times..."
	@previous_job_ids=""
	@i=1; \
	while [ $$i -le $(REPS) ]; do \
		echo "Iteration $$i"; \
		current_job_ids=""; \
		for job in $(SLURM_FILES); do \
			if [ -n "$$previous_job_ids" ]; then \
				jobid=`sbatch --dependency=afterok:$$previous_job_ids $$job | awk '{print $$4}'`; \
			else \
				jobid=`sbatch $$job | awk '{print $$4}'`; \
			fi; \
			current_job_ids="$$current_job_ids:$$jobid"; \
			echo "Submitted $$job with Job ID $$jobid"; \
		done; \
		previous_job_ids=$$(echo $$current_job_ids | sed 's/^://'); \
		i=$$((i+1)); \
	done

# Target to execute plot.slurm specifically
plot: $(PLOT_SLURM)
	@echo "Submitting plot SLURM job..."
	sbatch $(PLOT_SLURM)

# Target to clean only generated files
clean:
	@echo "Cleaning generated files..."
	rm -rf csv_aggregated csv_per_iterations logs plots

# Phony targets
.PHONY: all execute plot clean
