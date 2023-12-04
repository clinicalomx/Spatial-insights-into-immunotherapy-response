# usage > python Script_7_tTest.py --input Tumor.csv --output results.csv --group1 Responder --group2 Non_Responder --label Group

import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import argparse

def perform_t_tests(df, group1, group2, label_column, alpha=0.05):
    result_columns = ['Feature', 'Group1', 'Group2', 'T-Statistic', 'P-Value', 'Fold Change', 'Adjusted P-Value']
    result_df = pd.DataFrame(columns=result_columns)
    all_p_values = []  # Collect all p-values for later correction

    for feature in df.columns.difference([label_column, 'Image']):
        p_values = []

        group1_data = pd.to_numeric(df[df[label_column] == group1][feature], errors='coerce')
        group2_data = pd.to_numeric(df[df[label_column] == group2][feature], errors='coerce')

        t_statistic, p_value = ttest_ind(group1_data, group2_data, nan_policy='omit')
        p_values.append(p_value)

        fold_change = group1_data.mean() / group2_data.mean()

        result_df = result_df.append({
            'Feature': feature,
            'Group1': group1,
            'Group2': group2,
            'T-Statistic': t_statistic,
            'P-Value': p_value,
            'Fold Change': fold_change,
            'Adjusted P-Value': None
        }, ignore_index=True)

        all_p_values.extend(p_values)  # Collect p-values for correction

    print(f"Original P-values: {all_p_values}")

    reject, adjusted_p_values, _, _ = multipletests(all_p_values, alpha=alpha, method='fdr_bh')

    print(f"Adjusted P-values: {adjusted_p_values}")

    # Assign adjusted p-values to the corresponding rows in the result dataframe
    result_df['Adjusted P-Value'] = adjusted_p_values

    return result_df

def main():
    parser = argparse.ArgumentParser(description='Perform t-tests and save results to a CSV file.')
    parser.add_argument('--input', type=str, help='Input CSV file path')
    parser.add_argument('--output', type=str, help='Output CSV file path')
    parser.add_argument('--group1', type=str, help='Name of the first group')
    parser.add_argument('--group2', type=str, help='Name of the second group')
    parser.add_argument('--label', type=str, help='Column representing the group label')

    args = parser.parse_args()

    df = pd.read_csv(args.input)
    result_df = perform_t_tests(df, args.group1, args.group2, args.label)

    result_df.to_csv(args.output, index=False)

    print(result_df)

if __name__ == "__main__":
    main()
