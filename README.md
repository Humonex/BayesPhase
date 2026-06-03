# BayesPhase

BayesPhase is a research-stage Python tool for extending phased variant blocks using long-read methylation and SNP evidence. The current implementation bridges adjacent phased blocks from an input VCF and BAM file, with an optional joint methylation plus SNP phasing mode.

> **Pre-publication notice**
>
> This repository contains unpublished research code. Please keep the repository private until the associated manuscript is accepted or released as a preprint.

## Features

- Reads phased variants from a VCF file.
- Uses long-read alignments from a BAM file to evaluate evidence across phased blocks.
- Supports methylation-based block bridging.
- Supports optional joint methylation and SNP read phasing via `--jointPhase`.
- Writes an extended phased VCF as output.


## Requirements

Recommended environment:

- Python 3.9 or later
- `pysam`
- `pandas`
- `scipy`
- `cliffs_delta`

Install dependencies in a dedicated environment:

```bash
conda create -n bayesphase python=3.10
conda activate bayesphase
pip install pysam pandas scipy cliffs-delta
```

If `cliffs-delta` is not available in your package index, install the equivalent package used in your development environment that provides:

```python
from cliffs_delta import cliffs_delta
```

## Usage

```bash
python BayesPhase_joint_phase.py <input.vcf.gz> <input.bam> <output.vcf> \
  --threads 8 \
  --jointPhase
```

Arguments:

- `vcf`: Input phased VCF file.
- `bam`: Input BAM file containing long-read alignments.
- `out_vcf`: Output VCF path.
- `-t, --threads`: Number of worker processes. Default: `1`.
- `-jointPhase, --jointPhase`: Enable joint methylation plus SNP read phasing.

Example:

```bash
python BayesPhase_joint_phase.py \
  sample.phased.vcf.gz \
  sample.alignments.bam \
  sample.bayesphase.vcf \
  -t 8 \
  --jointPhase
```

## Input Files

The input VCF should contain phased genotype information and phase-set (`PS`) annotations. The input BAM should be indexed and should contain the read-level tags required by the workflow, including haplotype-related tags when available.

Example expected files:

```text
sample.phased.vcf.gz
sample.phased.vcf.gz.tbi
sample.alignments.bam
sample.alignments.bam.bai
```

## Output

BayesPhase writes an output VCF containing updated phase-set assignments after block bridging. The exact output interpretation should be described in the manuscript and validated against the benchmark or simulation protocol used in the study.

## Reproducibility Notes

Before public release, consider adding:

- A minimal test dataset.
- A small expected-output VCF.
- Exact software versions used in the manuscript.
- Benchmark scripts for reproducing reported figures and tables.
- A `requirements.txt` or `environment.yml`.

## Citation

Citation information will be added after the manuscript is available.

```bibtex
@article{bayesphase_tbd,
  title   = {BayesPhase: TODO},
  author  = {TODO},
  journal = {TODO},
  year    = {TODO}
}
```

## License

License information is not yet specified. Keep this repository private until the release policy is finalized.

## Contact

For questions about the method or implementation, contact:

- TODO: author name
- TODO: email address
