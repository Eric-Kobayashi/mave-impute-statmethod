# MaveDB Imputation Pipeline

This pipeline performs imputation on MaveDB single variant data using various simple statistical learning estimators.

Initial investigation leading to this pipeline is shown at [notebook](https://github.com/allydunham/xai_hackathon2024/blob/sklearn-imputation/mavedb-imputation.ipynb)  
Data analysis dashboard showing the pipeline results is hosted at a [temporary site](http://jira-report.hgi.sanger.ac.uk:8050/)

## Overview

The pipeline consists of two main scripts:
- `preprocess.py` - One-time data preprocessing
- `main.py` - Imputation for individual studies

## Quick Start

### 1. Preprocessing (one-time setup)
```bash
make preprocess
```

### 2. Submit batch jobs
```bash
make batch
```

## Detailed Usage

### Preprocessing
The preprocessing step:
- Reads raw `mavedb_singles_wide.tsv` data
- Removes URNs with high NaN values (>20%)
- Fills wild-type amino acid scores with zero
- Splits data by URN into individual files
- Generates `plan.txt` for batch submission

```bash
# Using default data directory
make preprocess

# Using custom data directory
make preprocess DATA_DIR=/path/to/data

# Direct script execution
python preprocess.py --data-dir ~/work/mavedb/data
```

### Batch Job Submission
The batch target:
- Runs preprocessing (if needed)
- Loads the wr module
- Submits all imputation jobs using `wr add -f plan.txt`

```bash
make batch
```

### Single Study Imputation
For testing or running individual studies:

```bash
# Random Forest
python main.py --urn urn:mavedb:00000001-a-1 --estimator-type random_forest --n-estimators 100

# KNN
python main.py --urn urn:mavedb:00000001-a-1 --estimator-type knn --n-neighbors 15

# Bayesian Ridge
python main.py --urn urn:mavedb:00000001-a-1 --estimator-type bayesian_ridge
```

## Supported Estimators

1. **Random Forest** (`random_forest`)
   - Parameter: `--n-estimators` (default: 100)

2. **K-Nearest Neighbors** (`knn`)
   - Parameter: `--n-neighbors` (default: 15)

3. **Bayesian Ridge** (`bayesian_ridge`)
   - No additional parameters

## Output Files

- **Preprocessed data**: `data/preprocessed/{urn}.tsv`
- **Imputed data**: `data/imputed/{urn}_{estimator}_{params}.tsv`
- **Evaluation results**: `data/evals/{urn}_{estimator}_{params}.tsv`

### Evaluation Output

The evaluation results contain performance metrics for each imputation model, calculated using a cross-validation approach:

- **MSE (Mean Squared Error)**: Measures the average squared difference between true and imputed values
- **Correlation**: Pearson correlation coefficient between true and imputed values
- **Spearman**: Spearman rank correlation coefficient between true and imputed values

Each evaluation file contains 100 rows, representing 100 iterations of the evaluation process where:
1. 10% of valid amino acid scores are randomly masked
2. The imputation model predicts these masked values
3. Performance metrics are calculated by comparing predictions to true values

This approach provides an estimate of model performance and accounts for the inherent variability in the imputation task.

## Makefile Targets

- `make help` - Show available targets and options
- `make preprocess` - Run data preprocessing
- `make batch` - Submit batch jobs using wr
- `make clean` - Remove generated files
- `make run-single` - Run single imputation job (for testing)

## Requirements

- Python 3.9+
- Required packages: pandas, numpy, scikit-learn, scipy, loguru, joblib, tqdm
    On the Sanger farm we have built an environment that contains all the dependencies. 
    /software/hgi/softpack/installs/users/eh19/mave-impute/1-scripts/python
- HGI wr module for batch submission 