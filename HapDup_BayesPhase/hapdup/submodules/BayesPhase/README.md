# BayesPhase Submodule

The local HapDup integration snapshot contained a `submodules/BayesPhase/` directory. In this repository, the maintained BayesPhase implementation is available at the repository root:

```text
BayesPhase_joint_phase.py
misc.py
```

The HapDup integration calls BayesPhase from `hapdup/main.py` after Margin phasing and WhatsHap retagging. Configure the `BayesPhase Path` placeholder in `hapdup/main.py` to point to the BayesPhase script used in your runtime environment.
