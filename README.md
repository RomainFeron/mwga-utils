# Collection of utilities for processing Multispecies Whole Genome Alignments

## Installation

Clone the repository and run make:

```bash
git clone https://github.com/RomainFeron/mgwa_utils.git  # Use ssh to clone if you have it setup
make -j 4  # -j specifies the number of threads to use in compilation
```

## Usage

### metrics

Compute basic metrics for each genomic position in the reference for a MAF file. The output is a wig file for each metric with a value for each genomic position in the reference assembly.

Available metrics:

- **Alignability**: proportion of assemblies aligned at this genomic position
- **Identity**: proportion of assemblies with the same nucleotide as the reference at this position

```
Usage:
    metrics <maf_file> [-p <prefix> -t <threads> -n <assemblies>]

Options:
    <maf_file>       Path to a MAF file.
    -p <prefix>      Prefix for output wig files [default: metrics]
    -n <assemblies>  Manually specify the number of assemblies in the alignment; if left to 0, it is computed from the MAF [default: 0]
    -t <threads>     Number of threads to use [default: 1].
    -h --help        Show this screen.
```

### missing_regions

Add regions from the reference assembly that are missing from a MAF file. Added regions have a score of **NA** and only contain the sequence from the reference assembly. The complete MAF file is output to `stdout`.

```
Usage:
    missing_regions <maf_file> <reference>

Options:
    <maf_file>       Path to a MAF file.
    <reference>      Path to a FASTA file for the reference assembly.
    -h --help        Show this screen.
```

### single_coverage

Check that the each sequence in the reference assembly is only covered once in a MAF. The output is a table of counts for each coverage value in each contig from the reference assembly.

```
Usage:
    single_coverage <maf_file> [-t <threads>]

Options:
    <maf_file>       Path to a MAF file.
    -h --help        Show this screen.
```

### stats

Compute a alignment statistics for a MAF file. Currently implemented statistics:

- Number of BP aligned in each assembly

```
Usage:
    stats <maf_file> [-p <prefix>]

Options:
    <maf_file>       Path to a MAF file.
    -p <prefix>      Prefix for output stats files [default: stats]
    -h --help        Show this screen.
```
