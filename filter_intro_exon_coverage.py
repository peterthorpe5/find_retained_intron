import pandas as pd
import argparse

def filter_intron_coverage(input_file, output_file="filtered_output.tsv"):
    """
    Filters intron data based on specific coverage criteria:
    - percent_intron_covered > 90%
    - upstream_exon_coverage > 90%
    - downstream_exon_coverage > 90%

    Keeps all columns including 'transcript_id'.

    Args:
        input_file (str): Path to the input TSV file.
        output_file (str): Path to save the filtered output.
    """
    # Load the dataset
    df = pd.read_csv(input_file, sep='\t')

    # Convert necessary columns to numeric
    df["percent_intron_covered"] = pd.to_numeric(df["percent_intron_covered"], errors='coerce')
    df["upstream_exon_coverage"] = pd.to_numeric(df["upstream_exon_coverage"], errors='coerce')
    df["downstream_exon_coverage"] = pd.to_numeric(df["downstream_exon_coverage"], errors='coerce')

    # Apply filtering conditions
    filtered_df = df[(df["percent_intron_covered"] > 90) &
                     (df["upstream_exon_coverage"] > 90) &
                     (df["downstream_exon_coverage"] > 90)]

    # Keep all columns including 'transcript_id' and remove duplicates
    filtered_df = filtered_df.drop_duplicates()

    # Save the final results
    filtered_df.to_csv(output_file, sep='\t', index=False)
    print(f"Filtered data saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter intron retention data based on coverage thresholds.")
    parser.add_argument("-i", "--input_file", required=True, help="Path to the input TSV file.")
    parser.add_argument("-o", "--output_file", default="filtered_output.tsv", help="Path to save the filtered results (default: filtered_output.tsv)")

    args = parser.parse_args()

    filter_intron_coverage(args.input_file, args.output_file)
