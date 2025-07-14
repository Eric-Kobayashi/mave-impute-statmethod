#! /software/hgi/softpack/installs/users/eh19/mave-impute/1-scripts/python
import pandas as pd
import pathlib
import argparse
import tqdm
from loguru import logger

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


def remove_high_nan_urn(single_wide_df, thresh=0.2):
    missing_vals = (
        single_wide_df.groupby(["score_urn"])[amino_acid_cols]
        .apply(lambda x: x.isna().mean())
        .mean(axis=1)
        .sort_values()
    )
    return single_wide_df[
        single_wide_df["score_urn"].isin(missing_vals[missing_vals < thresh].index)
    ]


def fill_wt_an_with_zero(df):
    # for each row, if norm_synonymous_score is nan, then fill the wt corresponding norm score to zero
    df2 = []
    for _, row in df.iterrows():
        new_row = row.copy()
        if pd.isna(row["norm_synonymous_score"]):
            new_row[row["wt"] + "_norm_score"] = 0
        df2.append(new_row)
    return pd.DataFrame(df2)


def preprocess(base_dir):
    """Preprocess the raw data and split by URN."""
    logger.info(f"Reading data from {base_dir / 'mavedb_singles_wide.tsv'}")
    df = pd.read_csv(base_dir / "mavedb_singles_wide.tsv", sep="\t")

    logger.info("Removing URNs with high NaN values")
    df = remove_high_nan_urn(df)

    logger.info("Filling wild-type amino acid scores with zero")
    df = fill_wt_an_with_zero(df)

    save_dir = base_dir / "preprocessed"
    save_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Saving preprocessed data for {df['score_urn'].nunique()} URNs")
    for urn, sample_df in tqdm.tqdm(df.groupby("score_urn")):
        sample_df.to_csv(save_dir / f"{urn}.tsv", sep="\t", index=False)

    logger.info("Generating imputation plan for wr submission")
    with open("plan.txt", "w") as f:
        for urn in df["score_urn"].unique():
            for n_estimators in [4, 10, 40, 100]:
                f.write(
                    f"/nfs/users/nfs_e/eh19/work/mavedb/main.py --urn {urn} --estimator-type random_forest --n-estimators {n_estimators}\n"
                )
            for n_neighbors in [5, 10, 15, 20]:
                f.write(
                    f"/nfs/users/nfs_e/eh19/work/mavedb/main.py --urn {urn} --estimator-type knn --n-neighbors {n_neighbors}\n"
                )
            f.write(
                f"/nfs/users/nfs_e/eh19/work/mavedb/main.py --urn {urn} --estimator-type bayesian_ridge\n"
            )

    logger.info("Preprocessing complete!")


def main():
    parser = argparse.ArgumentParser(
        description="Preprocess MaveDB data for imputation"
    )
    parser.add_argument(
        "--data-dir",
        type=str,
        default="~/work/mavedb/data",
        help="Directory containing the raw data (default: ~/work/mavedb/data)",
    )
    args = parser.parse_args()

    base_dir = pathlib.Path(args.data_dir).expanduser()
    preprocess(base_dir)


if __name__ == "__main__":
    main()
