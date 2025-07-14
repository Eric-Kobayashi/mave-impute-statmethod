# Makefile for MaveDB imputation pipeline

# Default data directory
DATA_DIR ?= ~/work/mavedb/data

# Python executable
PYTHON ?= /software/hgi/softpack/installs/users/eh19/mave-impute/1-scripts/python

# Default values for imputation parameters
URN ?= urn:mavedb:00000001-a-1
ESTIMATOR_TYPE ?= bayesian_ridge
N_ESTIMATORS ?= 100
N_NEIGHBORS ?= 15

.PHONY: help preprocess batch clean

help:
	@echo "Available targets:"
	@echo "  preprocess  - Run data preprocessing (one-time setup)"
	@echo "  batch       - Submit imputation jobs using wr"
	@echo "  clean       - Remove generated files"
	@echo ""
	@echo "Variables:"
	@echo "  URN             - Study URN to process (default: urn:mavedb:00000001-a-1)"
	@echo "  DATA_DIR        - Data directory (default: ~/work/mavedb/data)"
	@echo "  ESTIMATOR_TYPE  - Estimator type for single run (default: bayesian_ridge)"
	@echo "  N_ESTIMATORS    - Number of estimators for RF (default: 100)"
	@echo "  N_NEIGHBORS     - Number of neighbors for KNN (default: 15)"
	@echo ""
	@echo "Examples:"
	@echo "  make preprocess"
	@echo "  make batch"
	@echo "  make preprocess DATA_DIR=/path/to/data"
	@echo "  python main.py --urn urn:mavedb:00000001-a-1 --estimator-type knn --n-neighbors 10 --n-jobs -1"

preprocess:
	@echo "Running data preprocessing..."
	$(PYTHON) preprocess.py --data-dir $(DATA_DIR)
	@echo "Preprocessing complete!"

batch: 
	@echo "Checking if plan.txt exists..."
	@if [ ! -f plan.txt ]; then \
		echo "Error: plan.txt not found. Please run 'make preprocess' first to generate the job plan."; \
		exit 1; \
	fi
	@echo "Submitting imputation jobs using wr..."
	/software/hgi/installs/wr/latest/wr add -f plan.txt
	@echo "Jobs submitted successfully!"

clean:
	@echo "Cleaning generated files..."
	rm -f plan.txt
	rm -rf $(DATA_DIR)/preprocessed
	rm -rf $(DATA_DIR)/imputed
	rm -rf $(DATA_DIR)/evals
	@echo "Clean complete!"

# Example single run target (for testing)
run-single:
	@echo "Running single imputation job..."
	$(PYTHON) main.py --urn $(URN) --estimator-type $(ESTIMATOR_TYPE) --n-estimators $(N_ESTIMATORS) --n-neighbors $(N_NEIGHBORS) --n-jobs -1 --data-dir $(DATA_DIR) 