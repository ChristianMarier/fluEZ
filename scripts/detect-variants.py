import pandas as pd
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='Detect low-frequency variants from input TSV.')
    parser.add_argument('--input_mp', type=str, required=True, help='Path to the input TSV file')
    parser.add_argument('--output_hf', type=str, required=True, help='Path to the output file for high-frequency variants')
    parser.add_argument('--output_lf', type=str, required=True, help='Path to the output file for low-frequency variants')
    return parser.parse_args()

def main():
    args = parse_arguments()

    # Read the input TXT file into a pandas dataframe
    df = pd.read_csv(args.input_mp, sep='\t')

    # Remove characters preceding and including the ":" symbol in the "VALUES" column
    df['VALUES'] = df['VALUES'].apply(lambda x: x.split(':', 1)[1] if ':' in x else x)
    df['VALUES'] = df['VALUES'].apply(lambda x: x.split(':', 1)[1] if ':' in x else x)
    df['VALUES'] = df['VALUES'].apply(lambda x: x.split(':', 1)[1] if ':' in x else x)
    df['VALUES'] = df['VALUES'].apply(lambda x: x.split(':', 1)[1] if ':' in x else x)
    df['VALUES'] = df['VALUES'].apply(lambda x: x.split(':', 1)[1] if ':' in x else x)

    # Add new columns for A1, A2, A3, and A4 (add extra A5 for rare example of 5th allele/ambiguous assignment)
    df[['A1', 'A2', 'A3', 'A4']] = df['VALUES'].str.split(',', expand=True, n=3)

    # Replace NaN with 0 in A3 and A4
    df['A3'] = df['A3'].fillna(0)
    df['A4'] = df['A4'].fillna(0)
    #df['A5'] = df['A5'].fillna(0)

    # Convert the columns to numeric
    df['A1'] = pd.to_numeric(df['A1'])
    df['A2'] = pd.to_numeric(df['A2'])
    df['A3'] = pd.to_numeric(df['A3'])
    df['A4'] = pd.to_numeric(df['A4'])
    #df['A5'] = pd.to_numeric(df['A5'])

    # Fix NA segment name from getting converted to NA value
    df['CHROM'] = df['CHROM'].fillna("NA") 

    # Create the SUM column
    df['SUM'] = df['A1'] + df['A2'] + df['A3'] + df['A4']

    # Remove rows where SUM equals 0
    df = df[df['SUM'] != 0]

    # Add columns for allele frequencies
    df['AF1'] = (df['A1'] / df['SUM']) * 100
    df['AF2'] = (df['A2'] / df['SUM']) * 100
    df['AF3'] = (df['A3'] / df['SUM']) * 100
    df['AF4'] = (df['A4'] / df['SUM']) * 100

    # Add new columns for highest frequency and classifications
    df['HIGHEST FREQ'] = df[['AF1', 'AF2', 'AF3', 'AF4']].max(axis=1)
    df['HF CLASS'] = df.apply(lambda row: 'REF' if row['HIGHEST FREQ'] == row['AF1'] else 'ALT', axis=1)

    # Add new columns for 2nd highest frequency and classifications
    df['2ND HF'] = df.apply(lambda row: sorted([row['AF1'], row['AF2'], row['AF3'], row['AF4']], reverse=True)[1], axis=1)
    df['2HF CLASS'] = df.apply(lambda row: 'REF' if row['2ND HF'] == row['AF1'] else 'ALT', axis=1)

    # Add new columns for 3rd highest frequency and classifications
    df['3RD HF'] = df.apply(lambda row: sorted([row['AF1'], row['AF2'], row['AF3'], row['AF4']], reverse=True)[2], axis=1)
    df['3HF CLASS'] = df.apply(lambda row: 'REF' if row['3RD HF'] == row['AF1'] else 'ALT', axis=1)

    # Add new columns for 4th highest frequency and classifications
    df['4TH HF'] = df.apply(lambda row: sorted([row['AF1'], row['AF2'], row['AF3'], row['AF4']], reverse=True)[3], axis=1)
    df['4HF CLASS'] = df.apply(lambda row: 'REF' if row['4TH HF'] == row['AF1'] else 'ALT', axis=1)

    # Add FromEnd column
    df['FromEnd'] = df['LENGTH'] - df['POS']

    # Filter out rows where POS <= 10 and FromEnd <= 10
    df = df[(df['POS'] > 10) & (df['FromEnd'] > 10)]

    # Extract high-frequency variants
    high_freq_df = df[df['HF CLASS'] == 'ALT'][['CHROM', 'POS', 'REF', 'ALT', 'A1', 'A2', 'A3', 'A4', 'SUM', 'AF1', 'AF2', 'AF3', 'AF4', 'HIGHEST FREQ', 'HF CLASS', '2ND HF', '2HF CLASS', '3RD HF', '3HF CLASS', '4TH HF', '4HF CLASS']]
    high_freq_df.to_csv(args.output_hf, sep='\t', index=False)

    # Extract low-frequency variants
    low_freq_df = df[(df['2HF CLASS'] == 'ALT') & (df['2ND HF'] >= 2.5) & (df['SUM'] >= 400)][['CHROM', 'POS', 'REF', 'ALT', 'A1', 'A2', 'A3', 'A4', 'SUM', 'AF1', 'AF2', 'AF3', 'AF4', 'HIGHEST FREQ', 'HF CLASS', '2ND HF', '2HF CLASS', '3RD HF', '3HF CLASS', '4TH HF', '4HF CLASS']]
    low_freq_df = low_freq_df.append(df[(df['3HF CLASS'] == 'ALT') & (df['3RD HF'] >= 2.5) & (df['SUM'] >= 400)][['CHROM', 'POS', 'REF', 'ALT', 'A1', 'A2', 'A3', 'A4', 'SUM', 'AF1', 'AF2', 'AF3', 'AF4', 'HIGHEST FREQ', 'HF CLASS', '2ND HF', '2HF CLASS', '3RD HF', '3HF CLASS', '4TH HF', '4HF CLASS']])
    low_freq_df = low_freq_df.append(df[(df['4HF CLASS'] == 'ALT') & (df['4TH HF'] >= 2.5) & (df['SUM'] >= 400)][['CHROM', 'POS', 'REF', 'ALT', 'A1', 'A2', 'A3', 'A4', 'SUM', 'AF1', 'AF2', 'AF3', 'AF4', 'HIGHEST FREQ', 'HF CLASS', '2ND HF', '2HF CLASS', '3RD HF', '3HF CLASS', '4TH HF', '4HF CLASS']])
    low_freq_df.to_csv(args.output_lf, sep='\t', index=False)

if __name__ == "__main__":
    main()




