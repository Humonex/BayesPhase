# Integration Code Notes

This directory records the BayesPhase-specific changes needed to integrate BayesPhase into HapDup.

The full project snapshot is stored in:

```text
HapDup_BayesPhase/hapdup/
```

The implementation is based on the original HapDup project:

- https://github.com/KolmogorovLab/hapdup

HapDup depends on several third-party submodules and external binaries, including Flye, Margin, PEPPER, minimap2, samtools, and WhatsHap. Those components should be installed from their upstream sources.

## Main Modified Stage

The primary BayesPhase-specific HapDup modification is in `HapDup_BayesPhase/hapdup/hapdup/main.py` around the Margin phasing step:

1. Margin produces a phased VCF and a haplotagged BAM.
2. WhatsHap retags reads using the Margin phased VCF.
3. BayesPhase bridges phase blocks using the phased VCF and retagged BAM.
4. HapDup uses the BayesPhase bridge BAM for polishing and structural polishing.

See `hapdup_bayesphase_integration.patch` for a compact patch-style summary of the same integration point.

## Runtime Paths to Configure

The integrated workflow contains runtime placeholders that should be configured before running on a new system:

```python
SINGULARITY = "The path of SINGULARITY"
MARGIN_MIRROR = "The path of MARGIN/.sif"
PEPPER_MIRROR = "The path of PEPPER/.sif"
MARGIN_CONFIG_DIR = "path/hapdup/submodules/margin/params/phase"
```

The BayesPhase command in the integration is represented as:

```python
[PYTHON, "BayesPhase Path", margin_vcf, whatshap_haplotagged_sort_bam, bridge_vcf, bridge_haplotagged_bam, ...]
```

Replace `"BayesPhase Path"` with the BayesPhase script path used in your runtime environment, for example `submodules/BayesPhase/BayesPhase.py` or the repository-root `BayesPhase_joint_phase.py`.
