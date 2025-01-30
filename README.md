# Retained Intron Detection Script

This script identifies **retained introns** in long-read RNA sequencing data by:
1. **Parsing** a GTF file to extract **exons** and **introns**.
2. **Building** coordinate sets for each gene/transcript.
3. **Reading** a BAM file, extracting alignment blocks for each read.
4. **Checking** for introns that have ≥50% coverage by each read’s alignment blocks.
5. **Reporting** additional details, including coverage of flanking exons.

---

## Features
- **Automatic intron extraction** from a GTF file (requires intron features or explicit exons).
- **Accurate retained intron detection** via read alignment blocks, preventing false positives.
- **Flanking exon coverage** (upstream and downstream exon coverage %).
- **Automatic indexing** of introns (1-based) and exons (1-based).
- **PEP8-compliant** docstrings and structured logging (`debug.log`).

---

## Why Use This Script?
When analysing long-read RNA-seq data, one of the key goals is to **distinguish spliced reads** (introns removed) from **partially- or unspliced reads** (introns retained). Naive methods can incorrectly mark introns as retained just because a read’s start/end happens to span intron coordinates. **This script** solves that by considering the **actual alignment blocks** and ensuring that **≥50%** of the intron’s length is covered by contiguous mapping. If the read genuinely “skips” the intron (spliced out), it doesn’t register as retained. Retained introns are interesting to some people ....

---

## Requirements
- **Python 3.8+**  
- **pysam** for BAM parsing: `pip install pysam`
- **A GTF file** (with explicit **exon** and **intron** entries or a file from which introns can be inferred).
- **A sorted, indexed BAM file** aligned to the same reference used to produce the GTF.

---

## Preparing the Input Data
1. **Reference Genome & GTF**  
   - Ensure your **GTF** includes **`intron`** features or can be processed to add introns. For example, 
   if you only have exons, you can generate introns using [GenomeTools](http://genometools.org/) or `gffread`. Make sure you also 
   use the command to retain the ID. 
2. **Index the BAM**  
   - Use `samtools index yourfile.bam` to create an index (`.bai`).
   - map with minimap:
   # Map with minimap2 and sort the BAM file using 16 threads
    `minimap2 -t "$threads" -ax map-ont "$reference_genome" "${input_folder}/input.fastq" | \
    samtools view -@ "$threads" -bS | \
    samtools sort -@ "$threads" -o "$output_bam"`

    # Index the sorted BAM file using 16 threads
    `samtools index -@ "$threads" "$output_bam"`
3. **Check Chromosome Naming**  
   - If your GTF uses `Chr1` but your BAM uses `1`, consider renaming. The script does a  replacement of `Chr`, but ensure consistency.

---

## How to Run
```bash
python find_retained_introns.py \\
    --gtf /path/to/your_intron_annotated.gtf \\
    --bam /path/to/your_aligned_reads.bam \\
    --output retained_introns.tsv




## Output Explanation
The script creates a tab-delimited (.tsv) file with one row per retained intron event. Columns include:

read_id: Name/ID of the read in the BAM.
gene_id: Gene identifier from the GTF.
transcript_id: Transcript identifier from the GTF.
intron_number: 1-based index of the intron in the transcript.
intron_start, intron_end: Genomic coordinates of the intron.
read_start, read_end: Minimal & maximal covered positions of this read on the reference (based on alignment blocks).
read_position (if included): A combined chrom:start-end notation for read coverage.
read_length: The length of the read or query sequence.
overlap_length: How many positions in this intron are covered by the read.
percent_intron_covered: The fraction of the intron covered by the read (≥50% to be considered retained).
upstream_exon_coverage, downstream_exon_coverage: The percent coverage of exons flanking this intron (0–100).
upstream_exon_index, downstream_exon_index: 1-based indices of the flanking exons.



# test

A test bam can be found in the test/data folder. 
The arabidopsis gtf is not provided. 


gt gff3 -addintrons -retainids original.gtf > introns_defined.gtf
