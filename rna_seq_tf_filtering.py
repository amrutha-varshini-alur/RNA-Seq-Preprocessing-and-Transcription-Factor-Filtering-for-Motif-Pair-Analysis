import pandas as pd
import numpy as np
import requests
import time
import argparse

# Step 1: Load and clean the input file
def load_and_clean_file(file_path, output_path):
    try:
        # Attempt to read with tab delimiter
        df = pd.read_csv(file_path, sep='\t')
        print("Successfully read file with tab delimiter.")
    except:
        df = pd.read_csv(file_path, sep=',')
        print("Successfully read file with comma delimiter.")
    
    # Remove everything after the '.' in 'Gene_ID'
    df['gene_id'] = df['gene_id'].str.split('.').str[0]
    df.to_csv(output_path, sep='\t', index=False) 
    print(f"Cleaned file saved to {output_path}")
    return df

# Step 2: Map Ensembl IDs to gene symbols
def get_gene_symbol_ensembl(gene_id):
    url = f"http://rest.ensembl.org/lookup/id/{gene_id}?"
    headers = {"Content-Type": "application/json"}
    try:
        response = requests.get(url, headers=headers)
        if response.status_code == 200:
            return response.json().get('display_name')
        elif response.status_code == 400:
            print(f"Invalid Ensembl ID: {gene_id}")
        else:
            print(f"Failed request for Ensembl ID: {gene_id}, status code: {response.status_code}")
    except requests.RequestException as e:
        print(f"Request error for Ensembl ID {gene_id}: {e}")
    return None

def map_gene_symbols(df, output_path):
    df['gene_symbol'] = df['gene_id'].apply(lambda x: get_gene_symbol_ensembl(x) or "Unknown")
    time.sleep(0.1) 
    df.to_csv(output_path, index=False)
    print(f"Gene symbols mapped and saved to {output_path}")
    return df

# Step 3: Merge TPM datasets and calculate log fold change
def merge_and_filter_tpm(huvec_path, imr90_path, merged_output_path, sorted_output_path):
    huvec_df = pd.read_csv(huvec_path).rename(columns={'TPM': 'TPM_HUVEC'})
    imr90_df = pd.read_csv(imr90_path).rename(columns={'TPM': 'TPM_IMR90'})
    merged_df = pd.merge(huvec_df[['Gene name', 'TPM_HUVEC']], 
                         imr90_df[['Gene name', 'TPM_IMR90']], 
                         on='Gene name')
    filtered_df = merged_df[
        (merged_df['TPM_HUVEC'] > 0) & 
        (merged_df['TPM_IMR90'] > 0) &
        ((merged_df['TPM_HUVEC'] >= 1) | (merged_df['TPM_IMR90'] >= 1))
    ].copy()
    filtered_df['Abs_log_fold_change'] = np.abs(np.log(filtered_df['TPM_HUVEC'] / filtered_df['TPM_IMR90']))
    filtered_df.to_csv(merged_output_path, index=False)
    print(f"Filtered TPM dataset saved to {merged_output_path}")
    
    sorted_df = filtered_df.sort_values(by='Abs_log_fold_change', ascending=False)
    sorted_df.to_csv(sorted_output_path, index=False)
    print(f"Sorted TPM dataset saved to {sorted_output_path}")
    return sorted_df

# Step 4: Filter TF data for the top 75 genes
def filter_tf_data(tf_path, sorted_tpm_path, output_path):
    tf_data = pd.read_csv(tf_path, index_col=0)
    sorted_genes = pd.read_csv(sorted_tpm_path)
    top_75_genes = sorted_genes['Gene name'][:75].tolist()
    filtered_tf_data = tf_data[top_75_genes]
    filtered_tf_data.to_csv(output_path)
    print(f"Filtered TF data saved to {output_path}")

# Main execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="RNA-Seq Data Processing Pipeline")
    parser.add_argument("--input_file", required=True, help="Path to the input RNA-Seq file")
    parser.add_argument("--huvec_file", required=True, help="Path to the HUVEC dataset")
    parser.add_argument("--IMR90_file", required=True, help="Path to the IMR90 dataset")
    parser.add_argument("--tf_file", required=True, help="Path to the TF dataset")
    parser.add_argument("--cleaned_file", required=True, help="Path to save the cleaned RNA-Seq file")
    parser.add_argument("--symbol_mapped_file", required=True, help="Path to save the symbol-mapped file")
    parser.add_argument("--merged_tpm_file", required=True, help="Path to save the merged TPM file")
    parser.add_argument("--sorted_tpm_file", required=True, help="Path to save the sorted TPM file")
    parser.add_argument("--filtered_tf_file", required=True, help="Path to save the filtered TF data")

    args = parser.parse_args()

    # Step 1: Clean the input file
    cleaned_df = load_and_clean_file(args.input_file, args.cleaned_file)
    
    # Step 2: Map gene symbols
    symbol_mapped_df = map_gene_symbols(cleaned_df, args.symbol_mapped_file)
    
    # Step 3: Merge and filter TPM datasets
    sorted_tpm_df = merge_and_filter_tpm(
        args.huvec_file,
        args.symbol_mapped_file,  
        args.merged_tpm_file,
        args.sorted_tpm_file
    )
    
    # Step 4: Filter TF data
    filter_tf_data(
        args.tf_file,
        args.sorted_tpm_file,
        args.filtered_tf_file
    )
