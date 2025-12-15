# SAMREADER

`samreader.py` is a Python tool that parses a SAM mapping file and analyses its content in both quality and quantity. This work was carried out as part of a Master's degree project in Bioinformatics at the University of Montpellier.

## Features

- **Read Filtering**:
  Select reads based on:
  - Mapping quality threshold (MAPQ filtering)
  - Fully mapped status (based on FLAG and CIGAR string)

- **Chromosome-Based Analysis**:
  - Counts mapped and unmapped reads per chromosome
  - Computes the distribution of reads across reference sequences

- **Alignment Statistics**:
  - Computes alignment length statistics (minimum, maximum, mean)
  - Identifies short, intermediate, and long alignments
  - Calculates the percentage of reads containing at least one indel

- **Pair and Orientation Analysis**:
  - Estimates the percentage of properly paired reads
  - Estimates the percentage of properly oriented read pairs (FR or RF)

- **Window-Based Coverage Analysis**:
  - Divides reference sequences into fixed-size windows
  - Computes read coverage per window
  - Computes the mean MAPQ per window

- **Visualization**:
  - Generates coverage plots along each chromosome
  - Colors coverage by mean MAPQ values for intuitive quality assessment

- **Export Results**:
  - Saves a detailed summary table in a text file
  - Exports one coverage plot per chromosome as PNG images


## Usage

    python3 samreader.py path/to/file.sam

## Author

Copyright © 2025 -- Thomas JUILLAC 
<thomas.juillac@etu.umontpellier.fr> <thomas.juillac@supagro.fr>


## Requirements

Requirements:
- Python 3
- matplotlib

Install required package:
```bash
pip install matplotlib
```
