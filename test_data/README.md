# Example Input

This directory is reserved for the small BayesPhase example input set.

The example input prepared for this project contains four files:

```text
test_snp.gz.vcf.gz
test_snp.gz.tbi
test_reads.bam
test_reads.bam.bai
```

File roles:

- `test_snp.gz.vcf.gz`: bgzip-compressed, tabix-indexed input VCF containing phased SNP records.
- `test_snp.gz.tbi`: tabix index for the input VCF.
- `test_reads.bam`: coordinate-sorted long-read BAM used by BayesPhase.
- `test_reads.bam.bai`: BAM index for `test_reads.bam`.

Run the example from the repository root:

```bash
python BayesPhase_joint_phase.py \
  test_data/test_snp.gz.vcf.gz \
  test_data/test_reads.bam \
  test_data/test_output.vcf \
  -t 4 \
  --jointPhase
```

Keep the index files next to their corresponding data files. The example is intended as a small command-line smoke test, not as the full experimental benchmark.

Full experimental results are archived at:

- https://zenodo.org/records/21018164
