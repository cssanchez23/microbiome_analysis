#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import pandas as pd
from mprod import table2tensor
from mprod.dimensionality_reduction import TCAM

def parse_args():
    parser = argparse.ArgumentParser(
        description='Run TCAM analysis, given data organized by Subject_ID (column 1) and Timepoint (column 2).'
    )
    parser.add_argument(
        '-o', '--output_dir', metavar='output_dir', type=str, required=True,
        help='Output file path, e.g. Outputs/TCAM/.'
    )
    parser.add_argument(
        '-t', '--data_type', metavar='data_type', type=str, required=True,
        help='Data type, e.g., taxa, enzymes, or modules. Should be one word. Underscores accepted.'
    )
    parser.add_argument(
        '-i', '--input_table', metavar='input_table', type=str, required=True,
        help=(
            'Input feature table/metadata in tsv format. Rows are samples. Columns are features. '
            'Subject_ID and Timepoint should be columns 1 and 2.'
        )
    )
    return parser.parse_args()

def main():
    args = parse_args()

    print(f"Running TCAM analysis on {args.input_table} and outputting results to {args.output_dir}")

    file_path = args.input_table
    data_table = pd.read_csv(file_path, index_col=[0, 1], sep="\t", dtype={'Subject_ID': str})

    data_table.rename(columns={k: f"Feature_{e+1}" for e, k in enumerate(data_table.columns)}, inplace=True)
    data_tensor, map1, map3 = table2tensor(data_table, missing_flag=True)

    tca = TCAM()
    tca_trans = tca.fit_transform(data_tensor)

    tca_loadings = tca.mode2_loadings  # Obtain TCAM loadings
    tca_var = tca.explained_variance_ratio_ * 100  # % explained variation per TCA factor

    tca_df = pd.DataFrame(tca_trans)
    tca_df.rename(index=dict(map(reversed, map1.items())), inplace=True)  # Use the inverse of map1 to denote each row with Subject ID

    # Save the TCAM transformed data and explained variance
    os.makedirs(args.output_dir, exist_ok=True)
    tca_df.to_csv(os.path.join(args.output_dir, f'{args.data_type}_tcam_transform.csv'), index=True)
    tca_var.tofile(os.path.join(args.output_dir, f'{args.data_type}_tcam_variance_explained.csv'), sep=',')

if __name__ == "__main__":
    main()
