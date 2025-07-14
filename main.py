#! /software/hgi/softpack/installs/users/eh19/mave-impute/1-scripts/python
import pandas as pd
import pathlib
import numpy as np
import warnings
import argparse
from dataclasses import dataclass
from scipy.stats import spearmanr

from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.linear_model import BayesianRidge
from loguru import logger
from joblib import Parallel, delayed

amino_acid_cols = [
    "A_norm_score",
    "C_norm_score",
    "D_norm_score",
    "E_norm_score",
    "F_norm_score",
    "G_norm_score",
    "H_norm_score",
    "I_norm_score",
    "K_norm_score",
    "L_norm_score",
    "M_norm_score",
    "N_norm_score",
    "P_norm_score",
    "Q_norm_score",
    "R_norm_score",
    "S_norm_score",
    "T_norm_score",
    "V_norm_score",
    "W_norm_score",
    "Y_norm_score",
    "stop_norm_score",
]


def split_train_test(df, test_size=0.2):
    randomised_df = df.sample(frac=1)
    return (
        randomised_df.iloc[: int(len(randomised_df) * (1 - test_size))],
        randomised_df.iloc[int(len(randomised_df) * (1 - test_size)) :],
    )


@dataclass
class Stats:
    mse: float
    corr: float
    spearman: float


def train_model(
    train,
    test,
    estimator_type="random_forest",
    n_estimators=100,
    n_neighbors=15,
    n_jobs=1,
):
    amino_acid_cols_filtered = [
        col
        for col in amino_acid_cols
        if col != "stop_norm_score" and not train[col].isna().all()
    ]
    logger.info(f"Skipping {set(amino_acid_cols) - set(amino_acid_cols_filtered)}")

    # Configure estimator based on type
    if estimator_type == "random_forest":
        random_forest_params = {
            "max_depth": 10,
            "bootstrap": True,
            "max_samples": 0.5,
            "n_jobs": n_jobs,
            "random_state": 17,
        }
        estimator = RandomForestRegressor(
            n_estimators=n_estimators, **random_forest_params
        )
    elif estimator_type == "knn":
        estimator = KNeighborsRegressor(n_neighbors=n_neighbors, n_jobs=n_jobs)
    elif estimator_type == "bayesian_ridge":
        estimator = BayesianRidge()
    else:
        raise ValueError(
            f"Unknown estimator type: {estimator_type}. Supported types: random_forest, knn, bayesian_ridge"
        )

    imputer = IterativeImputer(random_state=17, estimator=estimator)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = imputer.fit(train[amino_acid_cols_filtered])

    return model, amino_acid_cols_filtered


def run_imputation(df, model, amino_acid_cols_filtered):
    transformed_all = model.transform(df[amino_acid_cols_filtered])
    df2 = df.copy()
    df2[amino_acid_cols_filtered] = transformed_all
    return df2


def evaluate_imputation_model(model, test, amino_acid_cols):
    """
    Evaluate imputation model by artificially masking 10% of values and measuring reconstruction quality.

    Strategy:
    1. Extract amino acid data from test set
    2. Create random mask for 10% of valid values
    3. Apply mask and run imputation
    4. Calculate metrics on masked positions

    Returns Stats with MSE, correlation, and Spearman correlation.
    """
    # Step 1: Extract original amino acid data
    original_data = test[amino_acid_cols].copy()

    # Step 2: Create artificial NaN mask for 10% of valid values
    flat_original = original_data.values.flatten()
    valid_indices = np.where(~np.isnan(flat_original))[0]
    n_to_mask = int(len(valid_indices) * 0.1)
    mask_indices = np.random.choice(valid_indices, size=n_to_mask, replace=False)

    # Step 3: Apply mask to create test data with artificial NaNs
    flat_masked = flat_original.copy()
    flat_masked[mask_indices] = np.nan
    masked_amino_data = flat_masked.reshape(original_data.shape)

    masked_test = test.copy().assign(
        **{col: masked_amino_data[:, i] for i, col in enumerate(amino_acid_cols)}
    )

    # Step 4: Run imputation
    imputed_values = model.transform(masked_test[amino_acid_cols])

    # Step 5: Extract true and imputed values at masked positions for evaluation
    true_at_mask = flat_original[mask_indices]
    imputed_at_mask = imputed_values.flatten()[mask_indices]

    # Remove any remaining NaNs (safety check)
    valid_mask = ~np.isnan(true_at_mask) & ~np.isnan(imputed_at_mask)
    true_eval = true_at_mask[valid_mask]
    imputed_eval = imputed_at_mask[valid_mask]

    # Step 6: Calculate evaluation metrics
    return Stats(
        mse=np.mean(np.square(true_eval - imputed_eval)),
        corr=np.corrcoef(true_eval, imputed_eval)[0, 1],
        spearman=spearmanr(true_eval, imputed_eval).statistic,
    )


def main():
    parser = argparse.ArgumentParser(description="Run imputation for a single study")
    parser.add_argument("--urn", type=str, required=True, help="Study URN to process")
    parser.add_argument(
        "--estimator-type",
        type=str,
        default="random_forest",
        choices=["random_forest", "knn", "bayesian_ridge"],
        help="Type of estimator to use for imputation",
    )
    parser.add_argument(
        "--n-estimators",
        type=int,
        default=100,
        help="Number of estimators for RandomForest (default: 100)",
    )
    parser.add_argument(
        "--n-neighbors",
        type=int,
        default=15,
        help="Number of neighbors for KNN estimator (default: 15)",
    )
    parser.add_argument(
        "--data-dir",
        type=str,
        default="~/work/mavedb/data",
        help="Directory containing preprocessed data (default: ~/work/mavedb/data)",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=1,
        help="Number of jobs for parallel processing (default: 1 for batch jobs, use -1 for single runs)",
    )
    args = parser.parse_args()

    base_dir = pathlib.Path(args.data_dir).expanduser()
    data_folder = base_dir / "preprocessed"

    if not data_folder.exists():
        raise FileNotFoundError(
            f"Preprocessed data directory not found: {data_folder}. "
            "Please run preprocess.py first."
        )

    sample_file = data_folder / f"{args.urn}.tsv"
    if not sample_file.exists():
        raise FileNotFoundError(
            f"Study data not found: {sample_file}. "
            "Please ensure the URN is correct and preprocessing has been completed."
        )

    sample_df = pd.read_csv(sample_file, sep="\t")
    logger.info(
        f"Imputing {args.urn} with {len(sample_df)} samples using {args.estimator_type} (n_jobs={args.n_jobs})"
    )

    try:
        train, test = split_train_test(sample_df)
        model, amino_acid_cols_filtered = train_model(
            train,
            test,
            args.estimator_type,
            args.n_estimators,
            args.n_neighbors,
            args.n_jobs,
        )
        imputed_df = run_imputation(sample_df, model, amino_acid_cols_filtered)

        # Save imputed data
        save_dir = base_dir / "imputed"
        save_dir.mkdir(parents=True, exist_ok=True)

        # Create filename based on estimator type and parameters
        if args.estimator_type == "random_forest":
            filename = f"{args.urn}_{args.estimator_type}_{args.n_estimators}.tsv"
        elif args.estimator_type == "knn":
            filename = f"{args.urn}_{args.estimator_type}_{args.n_neighbors}.tsv"
        else:  # bayesian_ridge
            filename = f"{args.urn}_{args.estimator_type}.tsv"

        imputed_df.to_csv(save_dir / filename, sep="\t", index=False)
        logger.info(f"Saved imputed data to {save_dir / filename}")

        # Save evaluation results
        save_dir = base_dir / "evals"
        save_dir.mkdir(parents=True, exist_ok=True)

        logger.info("Running evaluation...")
        all_evals = Parallel(n_jobs=args.n_jobs)(
            delayed(evaluate_imputation_model)(model, test, amino_acid_cols_filtered)
            for _ in range(100)
        )
        all_evals_df = pd.DataFrame(
            [
                {"mse": eval.mse, "corr": eval.corr, "spearman": eval.spearman}
                for eval in all_evals
            ]
        )
        all_evals_df.to_csv(save_dir / filename, sep="\t", index=False)
        logger.info(f"Saved evaluation results to {save_dir / filename}")

    except Exception as e:
        logger.exception(f"Error imputing {args.urn}: {e}")
        raise


if __name__ == "__main__":
    main()
