import pysam
import csv
import argparse
import re
import logging
from collections import defaultdict

# Configure logging
logging.basicConfig(
    filename="debug.log",
    level=logging.DEBUG,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

def expand_coordinates(ranges):
    """
    Expand coordinate ranges into a full set of positions for easy lookup.

    Args:
        ranges (list of tuple): List of (start, end) coordinates.

    Returns:
        set: Set of all positions within the provided ranges.
    """
    expanded = set()
    for start, end in ranges:
        expanded.update(range(start, end + 1))
    return expanded


def parse_gtf(gtf_file):
    """
    Parse a GTF file and extract exon and intron features grouped by transcript and gene.
    Expands exons and introns into full coordinate sets for easy matching.

    Args:
        gtf_file (str): Path to the GTF file.

    Returns:
        dict: A dictionary where keys are gene IDs, values contain exon and intron positions.
              Each entry has:
                'chrom': Chromosome name
                'transcripts': {
                    transcript_id: {
                        'exons': [(start, end), ...],
                        'introns': [(start, end), ...],
                        'exon_coords': set(...),
                        'intron_coords': set(...),
                        'exon_list': [(start, end), ...]
                    }
                }
    """
    gene_structure = defaultdict(
        lambda: {
            "chrom": "",
            "transcripts": defaultdict(
                lambda: {
                    "exons": [],
                    "introns": [],
                    "exon_coords": set(),
                    "intron_coords": set(),
                    "exon_list": []
                }
            )
        }
    )

    transcript_to_gene = {}

    print(f"Parsing GTF file: {gtf_file}")
    logging.info("Parsing GTF file: %s", gtf_file)

    with open(gtf_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                logging.warning("Skipping malformed GTF line: %s", line.strip())
                continue

            feature_type = fields[2]
            chrom = fields[0].replace("Chr", "")  # Normalise chromosome names
            start = int(fields[3])
            end = int(fields[4])
            attributes = fields[8]

            gene_id = None
            transcript_id = None

            for attr in attributes.split(";"):
                attr = attr.strip()
                if "geneID=" in attr:
                    gene_id = attr.split("=")[1]
                elif "ID=" in attr and feature_type == "mRNA":
                    transcript_id = attr.split("=")[1]
                elif "Parent=" in attr:
                    transcript_id = attr.split("=")[1]

            # Map transcripts to genes
            if gene_id and transcript_id:
                transcript_to_gene[transcript_id] = gene_id
                gene_structure[gene_id]["chrom"] = chrom

            # Assign features to the correct gene
            assigned_gene = transcript_to_gene.get(transcript_id, None)
            if not assigned_gene:
                logging.warning("Skipping feature with missing gene mapping: %s", line.strip())
                continue

            # Collect exons/introns
            if feature_type == "exon":
                gene_structure[assigned_gene]["transcripts"][transcript_id]["exons"].append((start, end))
                gene_structure[assigned_gene]["transcripts"][transcript_id]["exon_list"].append((start, end))
            elif feature_type == "intron":
                gene_structure[assigned_gene]["transcripts"][transcript_id]["introns"].append((start, end))

    # Expand exon/intron coords
    for gene, data in gene_structure.items():
        for transcript, structure in data["transcripts"].items():
            structure["exon_coords"] = expand_coordinates(structure["exons"])
            structure["intron_coords"] = expand_coordinates(structure["introns"])

    print(f"Finished parsing GTF file. Total genes processed: {len(gene_structure)}")
    logging.info("Finished parsing GTF file. Total genes processed: %d", len(gene_structure))

    return gene_structure


def coverage_percentage(read_start, read_end, region_start, region_end):
    """
    Calculate what percent of [region_start, region_end] is covered by [read_start, read_end].

    Args:
        read_start (int): The read's start coordinate (0-based).
        read_end (int): The read's end coordinate (0-based).
        region_start (int): Start of the region.
        region_end (int): End of the region.

    Returns:
        float: Percent coverage of the region by the read, rounded to one decimal.
    """
    overlap_start = max(read_start, region_start)
    overlap_end = min(read_end, region_end)
    overlap_len = max(0, overlap_end - overlap_start + 1)
    region_len = region_end - region_start + 1
    if region_len == 0:
        return 0.0
    return round((overlap_len / region_len) * 100, 1)


def find_flanking_exons(intron_index, exon_list):
    """
    Identify the upstream and downstream exons based on a sorted exon list.

    The intron i is typically between exons i and i+1.

    Args:
        intron_index (int): Index of the intron in sorted introns for a transcript.
        exon_list (list of tuple): Sorted list of (start, end) exons.

    Returns:
        (tuple or None, tuple or None): (upstream_exon, downstream_exon)
    """
    upstream_exon = None
    downstream_exon = None

    if intron_index < len(exon_list):
        upstream_exon = exon_list[intron_index]
    if intron_index + 1 < len(exon_list):
        downstream_exon = exon_list[intron_index + 1]

    return upstream_exon, downstream_exon


def blocks_to_set(read):
    """
    Convert read alignment blocks (CIGAR-based contiguous regions) to a set of covered positions.

    Args:
        read (pysam.AlignedSegment): A read from pysam.

    Returns:
        set: All reference positions the read covers (excludes spliced-out / 'N' segments).
    """
    # read.get_blocks() -> list of (ref_start, ref_end) intervals
    covered = set()
    blocks = read.get_blocks()
    for ref_start, ref_end in blocks:
        covered.update(range(ref_start, ref_end))  # ref_end is exclusive
    return covered



def find_retained_introns(bam_file, gene_structure, output_file):
    """
    Identify retained introns by checking read alignment blocks vs intron coordinates.
    Require >= 50% coverage to call intron retained.
    Calculate upstream/downstream exon coverage as well.
    Includes exon indices and a combined 'read_position' column.

    Args:
        bam_file (str): Path to the BAM file.
        gene_structure (dict): Parsed gene structure from GTF file.
        output_file (str): Path to the output TSV file.
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    retained_introns = []

    print(f"Processing BAM file: {bam_file}")
    logging.info("Processing BAM file: %s", bam_file)

    for read in bam.fetch():
        if read.is_unmapped:
            continue
        if read.cigartuples is None:
            continue

        # Convert alignment blocks to a set of covered positions
        covered_positions = blocks_to_set(read)
        if not covered_positions:
            continue

        chrom = bam.get_reference_name(read.reference_id).replace("Chr", "")
        read_id = read.query_name
        read_length = read.query_length

        read_start = min(covered_positions)
        read_end = max(covered_positions)
        read_position = f"{chrom}:{read_start}-{read_end}"  # new column

        logging.debug(
            "Read %s on %s covers %d positions from %d to %d",
            read_id, chrom, len(covered_positions), read_start, read_end
        )

        # For each gene on this chromosome
        for gene_id, data in gene_structure.items():
            if data["chrom"] != chrom:
                continue

            # For each transcript of this gene
            for transcript_id, structure in data["transcripts"].items():
                sorted_introns = sorted(structure["introns"], key=lambda x: x[0])
                sorted_exons = sorted(structure["exon_list"], key=lambda x: x[0])

                for i, (i_start, i_end) in enumerate(sorted_introns):
                    intron_positions = range(i_start, i_end + 1)
                    intron_len = i_end - i_start + 1
                    # overlap with covered_positions
                    overlap = covered_positions.intersection(intron_positions)
                    overlap_len = len(overlap)

                    if overlap_len >= 0.5 * intron_len:  # require >=50% coverage
                        percent_intron_covered = round((overlap_len / intron_len) * 100, 1)

                        # flanking exons
                        upstream_exon, downstream_exon = find_flanking_exons(i, sorted_exons)
                        up_exon_cov = 0.0
                        upstream_exon_index = None
                        if upstream_exon:
                            # coverage
                            up_exon_cov = coverage_percentage(
                                read_start, read_end,
                                upstream_exon[0], upstream_exon[1]
                            )
                            # find which exon index in sorted_exons
                            for exon_idx, exon in enumerate(sorted_exons):
                                if exon == upstream_exon:
                                    upstream_exon_index = exon_idx + 1  # 1-based index
                                    break

                        down_exon_cov = 0.0
                        downstream_exon_index = None
                        if downstream_exon:
                            down_exon_cov = coverage_percentage(
                                read_start, read_end,
                                downstream_exon[0], downstream_exon[1]
                            )
                            for exon_idx, exon in enumerate(sorted_exons):
                                if exon == downstream_exon:
                                    downstream_exon_index = exon_idx + 1
                                    break

                        intron_number = i + 1

                        retained_introns.append({
                            "read_id": read_id,
                            "gene_id": gene_id,
                            "transcript_id": transcript_id,
                            "intron_number": intron_number,
                            "intron_start": i_start,
                            "intron_end": i_end,
                            "read_start": read_start,
                            "read_end": read_end,
                            "read_position": read_position,  # new column
                            "read_length": read_length,
                            "overlap_length": overlap_len,
                            "percent_intron_covered": percent_intron_covered,
                            "upstream_exon_coverage": up_exon_cov,
                            "downstream_exon_coverage": down_exon_cov,
                            "upstream_exon_index": upstream_exon_index,    # new column
                            "downstream_exon_index": downstream_exon_index # new column
                        })

    # Write output
    fieldnames = [
        "read_id", "gene_id", "transcript_id", "intron_number",
        "intron_start", "intron_end",
        "read_start", "read_end", "read_position",
        "read_length", "overlap_length", "percent_intron_covered",
        "upstream_exon_coverage", "downstream_exon_coverage",
        "upstream_exon_index", "downstream_exon_index"
    ]
    with open(output_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(retained_introns)

    print(f"Retained introns written to {output_file}")
    logging.info("Retained introns written to %s", output_file)


def main():
    """
    Main function to parse arguments and execute retained intron detection.
    """
    parser = argparse.ArgumentParser(
        description="Identify retained introns from BAM and GTF files."
    )
    parser.add_argument(
        "--gtf",
        type=str,
        default="genomic_data/Araport11_fixed_with_introns.gtf",
        help="Path to the GTF file with explicitly defined exons and introns."
    )
    parser.add_argument(
        "--bam",
        type=str,
        default="test.bam",
        help="Path to the BAM file. Default: test.bam"
    )
    parser.add_argument(
        "--output",
        type=str,
        default="retained_introns.tsv",
        help="Path to the output file. Default: retained_introns.tsv"
    )

    args = parser.parse_args()
    logging.info("Starting analysis")

    gene_structure = parse_gtf(args.gtf)
    find_retained_introns(args.bam, gene_structure, args.output)

    logging.info("Analysis complete. Results saved to %s", args.output)


if __name__ == "__main__":
    main()
