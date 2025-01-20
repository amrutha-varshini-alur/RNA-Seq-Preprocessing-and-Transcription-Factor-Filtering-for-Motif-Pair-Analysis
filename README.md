# RNA-Seq-Preprocessing-and-Transcription-Factor-Filtering-for-Motif-Pair-Analysis
This script, rna_seq_tf_filtering.py, is designed to facilitate preprocessing steps for motif pair analysis, focusing on RNA-Seq data, Transcription per million (TPM) filtering, and transcription factor (TF) selection. The script integrates multiple data preparation steps into a streamlined pipeline, making it efficient for large-scale bioinformatics workflows.

## Features

The pipeline includes the following steps:

1. **Loading and Cleaning Input Files**
* Handles tab- or comma-delimited RNA-Seq data files.
* Cleans gene IDs by removing version numbers (everything after '.').
* Saves the cleaned file for downstream processing.

2. **Mapping Ensembl IDs to Gene Symbols**
* Uses Ensembl's REST API to map Ensembl IDs to their respective gene symbols.
* Fallback mechanism to handle errors or unmapped IDs, assigning them as "Unknown."
* Outputs a file with mapped gene symbols.

3. **Merging and Filtering TPM Datasets**
* Merges TPM (Transcripts Per Million) data from two cell lines (e.g., HUVEC and IMR90).
* Filters out genes with non-expressive or low-expression values.
* Calculates the absolute log fold change between TPM values.
* Outputs both filtered and sorted datasets.

4. **Filtering TF Data**
* Filters a transcription factor dataset to include only the top 75 genes based on absolute log fold change.
* Saves the filtered TF dataset for motif analysis.

## Prerequisites

**Python Libraries**

Ensure the following Python libraries are installed:
* pandas
* numpy
* requests
* argparse
  
To install the required libraries, run:<br>
pip install pandas numpy requests argparse

## Data Requirements

* **Input RNA-Seq file:** Should include a column gene_id with Ensembl IDs and associated TPM values.
* **TPM Files:** Files for HUVEC and IMR90 datasets, containing Gene name and TPM columns.
* **TF Dataset:** Contains transcription factor data with genes as columns.

## Usage

Run the script using the command line:

python motif_pair_preprocessing.py --input_file <path_to_input_file> \ <br>
                                   --huvec_file <path_to_huvec_file> \ <br>
                                   --IMR90_file <path_to_imr90_file> \ <br>
                                   --tf_file <path_to_tf_file> \ <br>
                                   --cleaned_file <path_to_cleaned_file> \ <br>
                                   --symbol_mapped_file <path_to_symbol_mapped_file> \ <br>
                                   --merged_tpm_file <path_to_merged_tpm_file> \ <br>
                                   --sorted_tpm_file <path_to_sorted_tpm_file> \ <br>
                                   --filtered_tf_file <path_to_filtered_tf_file> <br>

**Command-Line Arguments**
| Argument	             |                                              Description |
|:---------------------  |  :-------------------------------------------------------|
| --input_file	         |                                Path to the RNA-Seq file. |
| --huvec_file	         |                           Path to the HUVEC TPM dataset. |
| --IMR90_file	         |                           Path to the IMR90 TPM dataset. |
| --tf_file	             |                Path to the transcription factor dataset. |
| --cleaned_file	       |                Output path for the cleaned RNA-Seq file. |
| --symbol_mapped_file   |       Output path for the file with mapped gene symbols. |
| --merged_tpm_file	     |                  Output path for the merged TPM dataset. |
| --sorted_tpm_file	     |                  Output path for the sorted TPM dataset. |
| --filtered_tf_file	   |  Output path for the filtered transcription factor data. |

## Workflow Details

**Step 1: Load and Clean RNA-Seq File**
* Removes the version suffix in gene_id (e.g., ENSG000001.1 → ENSG000001).
* Handles both tab- and comma-delimited files.

**Step 2: Map Gene Symbols**
* Queries Ensembl’s REST API for each Ensembl ID to fetch gene symbols.
* Handles API errors gracefully with fallback to "Unknown."

**Step 3: Merge and Filter TPM Data**
* Merges HUVEC and IMR90 TPM datasets by Gene name.
* Filters genes based on:
  * Positive TPM values in both datasets.
  * Minimum TPM threshold of 1 in either dataset.
* Computes the absolute log fold change (|log(TPM_HUVEC / TPM_IMR90)|).
* Outputs filtered and sorted datasets.

**Step 4: Filter TF Data**
* Extracts top 75 genes based on log fold change.
* Filters TF dataset to include only the selected genes.

## Outputs

The script generates the following files:
1. Cleaned RNA-Seq File: Processed version of the input file with cleaned gene_id.
2. Gene Symbol Mapped File: RNA-Seq data with an additional column for gene_symbol.
3. Merged TPM Dataset: Combined HUVEC and IMR90 TPM data.
4. Sorted TPM Dataset: TPM data sorted by absolute log fold change.
5. Filtered TF Dataset: TF data containing only the top 75 genes.
