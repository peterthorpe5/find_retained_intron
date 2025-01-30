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



Output Explanation
The script generates a tab-separated values (TSV) file, where each row represents an intron retention event detected from the BAM file. The columns include:

Column Name	Description
read_id	The read name (or ID) from the BAM file.
gene_id	The gene associated with this transcript (from the GTF).
transcript_id	The transcript ID assigned from the GTF.
intron_number	The intron index in this transcript (1-based).
intron_start	The genomic start coordinate of the intron.
intron_end	The genomic end coordinate of the intron.
read_start	The genomic start coordinate of the read.
read_end	The genomic end coordinate of the read.
read_position	The formatted position of the read in the form chrom:start-end.
read_length	The total length of the read.
overlap_length	The number of bases the read overlaps with the intron.
percent_intron_covered	The percentage of the intron covered by the read (≥50% required to be considered retained).
upstream_exon_coverage	The percentage of the upstream exon covered by the read.
downstream_exon_coverage	The percentage of the downstream exon covered by the read.
upstream_exon_index	The 1-based index of the upstream exon in the transcript.
downstream_exon_index	The 1-based index of the downstream exon in the transcript.




# test

A test bam can be found in the test/data folder. 
The arabidopsis gtf is not provided. 


gt gff3 -addintrons -retainids original.gtf > introns_defined.gtf

# generate test: test with known "AT1G01010.1_unspliced" unspliced
minimap2 -t 4 -ax map-ont ./genomic_data/TAIR10_chr_all.fas spliced_unsplied.fa |     samtools view -@ "$threads" -bS |      samtools sort -@ "$threads" -o test.bam

samtools index test.bam


the results folder for the tests:


Below is a short summary of how these test read IDs relate to the intron retention outcomes you see in the table:

AT1G01010.1_UNSPLICED

This read is designed to retain all introns for the AT1G01010 gene.
Consequently, it shows 100% coverage of every intron (intron numbers 1 to 5).
AT1G01010.1_INTRON_1_2_WITH_EXONS/19-975

This read is configured to have introns 1 and 2 retained while leaving other introns spliced out.
As a result, the table shows high coverage for introns 1 and 2 but not for the subsequent ones.
AT1G01010.1_EXON1_2_INTRON1/22-642

This read contains exons 1 and 2 but retains the first intron only.
So the table lists coverage specifically for intron 1 (91.4–100%) while other introns aren’t retained.
Because these reads are contrived test cases, their names reflect the expected structure (splicing vs. retention). The retained intron detection in the output table confirms whether each test read matches its intended intron-retention profile.