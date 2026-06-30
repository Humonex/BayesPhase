# Integration Code Notes

This directory records the BayesPhase-specific changes needed to integrate BayesPhase into HapDup.

The local implementation was inspected from:

```text
C:/Users/luomo/Nutstore/1/我的坚果云/BayesPhase/code/HapDup_BayesPhase/hapdup
```

The implementation is intentionally stored here as a patch-style integration note instead of vendoring the full HapDup source tree. HapDup depends on several third-party submodules and external binaries, including Flye, Margin, PEPPER, minimap2, samtools, and WhatsHap. Those components should be installed from their upstream sources.

## Main Modified Stage

The integration modifies HapDup's `hapdup/main.py` around the Margin phasing step:

1. Margin produces a phased VCF and a haplotagged BAM.
2. WhatsHap retags reads using the Margin phased VCF.
3. BayesPhase bridges phase blocks using the phased VCF and retagged BAM.
4. HapDup uses the BayesPhase bridge BAM for polishing and structural polishing.

See `hapdup_bayesphase_integration.patch` for the implementation sketch.

## Runtime Paths to Configure

The local implementation uses placeholder paths that should be configured before running on a new system:

```python
SINGULARITY = "The path of SINGULARITY"
MARGIN_MIRROR = "The path of MARGIN/.sif"
PEPPER_MIRROR = "The path of PEPPER/.sif"
MARGIN_CONFIG_DIR = "path/hapdup/submodules/margin/params/phase"
```

The BayesPhase command in the local implementation is represented as:

```python
[PYTHON, "BayesPhase Path", margin_vcf, whatshap_haplotagged_sort_bam, bridge_vcf, bridge_haplotagged_bam, ...]
```

Replace `"BayesPhase Path"` with the actual BayesPhase script path used in your environment.
