# BayesPhase

BayesPhase is a research-stage Python tool for extending phased variant blocks using long-read methylation and SNP evidence. The current implementation bridges adjacent phased blocks from an input VCF and BAM file, with an optional joint methylation plus SNP phasing mode.


## Features

- Reads phased variants from a VCF file.
- Uses long-read alignments from a BAM file to evaluate evidence across phased blocks.
- Supports methylation-based block bridging.
- Supports optional joint methylation and SNP read phasing via `--jointPhase`.
- Writes an extended phased VCF as output.

## Installation

The recommended setup is the conda environment in `environment.yml`, because `pysam` depends on the HTSlib ecosystem used for BAM/VCF access.

```bash
conda env create -f environment.yml
conda activate bayesphase
```

A pip-only dependency list is also provided for environments where the required system libraries are already available:

```bash
python -m pip install -r requirements.txt
```

Core Python dependencies:

- Python 3.10 or later
- `pysam`
- `pandas`
- `scipy`
- `cliffs-delta` (`import cliffs_delta`)

## Usage

```bash
python BayesPhase.py <input.vcf.gz> <input.bam> <output.vcf> \
  --threads 8 \
  --jointPhase
```

Arguments:

- `vcf`: Input phased VCF file.
- `bam`: Input BAM file containing long-read alignments.
- `out_vcf`: Output VCF path.
- `-t, --threads`: Number of worker processes. Default: `1`.
- `-jointPhase, --jointPhase`: Enable joint methylation plus SNP read phasing.

Methylation-only mode can be run by omitting `--jointPhase`.

## Input Files

The input VCF should be bgzip-compressed, tabix-indexed, and contain phased genotype information with phase-set (`PS`) annotations. The input BAM should be coordinate-sorted, indexed, and contain the read-level tags required by the workflow, including methylation tags and haplotype-related tags when available.

Expected file pairs:

```text
sample.phased.vcf.gz
sample.phased.vcf.gz.tbi
sample.alignments.bam
sample.alignments.bam.bai
```

## Example Input

A small example input set is documented in `test_data/`. The files prepared for the example are:

```text
test_data/test_snp.gz.vcf.gz
test_data/test_snp.gz.tbi
test_data/test_reads.bam
test_data/test_reads.bam.bai
```

Example command from the repository root:

```bash
python BayesPhase.py \
  test_data/test_snp.gz.vcf.gz \
  test_data/test_reads.bam \
  test_data/test_output.vcf \
  -t 4 \
  --jointPhase
```

See `test_data/README.md` and `test_data/MANIFEST.tsv` for the example input manifest.

## Output

BayesPhase writes an output VCF containing updated phase-set assignments after block bridging. The output VCF can be compressed and indexed with standard HTSlib tools if downstream workflows require bgzip/tabix files.

## Experimental Results

The experimental results associated with this project are archived on Zenodo:

- [https://zenodo.org/records/21018164](https://zenodo.org/records/21018164)

## Repository Layout

```text
BayesPhase.py               Main command-line implementation
misc.py                     Shared logging and file helpers
environment.yml            Recommended conda environment
requirements.txt           Pip dependency list
test_data/                 Example input manifest and example data location
HapDup_BayesPhase/          HapDup-based BayesPhase integration project snapshot
```


## License

This project is distributed under the GNU General Public License v3.0. 
